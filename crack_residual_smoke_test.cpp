/*
 Local build-verification test for CrackResidual2D, issue #16.

 Motivating case: two rigid pieces sliding past each other along a known
 vertical line, built from real speckle texture (examples/2d_dic/oht_cfrp_0.bmp)
 rather than a synthetic pattern, with an exact-integer-pixel shift so the
 ground truth has no sub-pixel interpolation bias.

 The point of this test is specifically to show what NEITHER per-POI ZNCC nor
 sigma/beta can see: every POI in the grid sits far enough from the line that
 its own subset window never straddles it, so each one is individually a
 perfectly good correlation. CrackResidual2D should still light up in a band
 centered on the line (evaluated via POIs from both sides, since its search
 radius is larger than any one subset), and stay low away from it.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;
	int thread_number = 4;
	omp_set_num_threads(thread_number);

	//--- build the two-piece "crack" image pair from real speckle texture ---
	cout << "=== Building synthetic two-rigid-piece image pair ===" << endl;
	Image2D source_img("examples/2d_dic/oht_cfrp_0.bmp");

	int x_lo = 30, x_hi = 270, y_lo = 150, y_hi = 350;
	int width = x_hi - x_lo, height = y_hi - y_lo;
	int x0 = 120; //crack line, in the cropped image's own local coordinates (0..width-1), centered
	int shift = 5; //exact-integer rigid shift on the right-hand side

	Image2D ref_img(width, height);
	Image2D tar_img(width, height);
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int global_x = x + x_lo, global_y = y + y_lo;
			ref_img.cv_mat.at<uchar>(y, x) = source_img.cv_mat.at<uchar>(global_y, global_x);

			int u = (x >= x0) ? shift : 0;
			int source_x = global_x - u; //target(x) = reference(x - u), matching ICGN's own warp convention
			tar_img.cv_mat.at<uchar>(y, x) = source_img.cv_mat.at<uchar>(global_y, source_x);
		}
	}
	cv::cv2eigen(ref_img.cv_mat, ref_img.eg_mat);
	cv::cv2eigen(tar_img.cv_mat, tar_img.eg_mat);

	//--- POI grid: columns far enough from x0 that no subset window straddles it ---
	int subset_radius = 8;
	vector<int> col_offsets = { -60, -40, -20, -12, 12, 20, 40, 60 };
	vector<POI2D> poi_queue;
	for (int y = 40; y <= height - 40; y += 20)
	{
		for (int offset : col_offsets)
		{
			poi_queue.push_back(POI2D(Point2D((float)(x0 + offset), (float)y)));
		}
	}

	//seed with a coarse FFTCC search before ICGN -- with a real 5px whole-pixel shift and a
	//zero initial guess, ICGN's gradient-based Newton iteration doesn't reliably converge to
	//the true answer on this real speckle texture (verified: it diverges to nonsense values on
	//most right-side POIs without this), the same reason other smoke tests in this repo
	//(e.g. simplex_smoke_test.cpp) always seed ICGN from FFTCC2D rather than a zero guess
	FFTCC2D fftcc(subset_radius, subset_radius, thread_number);
	fftcc.setImages(ref_img, tar_img);
	fftcc.compute(poi_queue);

	ICGN2D1 icgn(subset_radius, subset_radius, 0.001f, 20, thread_number);
	icgn.setImages(ref_img, tar_img);
	icgn.prepare();
	icgn.compute(poi_queue);

	int good_count = 0;
	for (auto& poi : poi_queue)
	{
		if (poi.result.zncc > 0.9f) good_count++;
	}
	cout << "  " << good_count << "/" << poi_queue.size() << " POIs individually well-correlated (zncc > 0.9)" << endl;
	bool all_good_ok = good_count == (int)poi_queue.size();
	cout << "  " << (all_good_ok ? "PASS" : "FAIL") << ": every POI's own subset avoided the line, so none show any sign of the discontinuity" << endl;
	if (!all_good_ok) failures++;

	//--- CrackResidual2D: should see what ZNCC can't ---
	cout << endl << "=== CrackResidual2D residual, near the line vs. far from it ===" << endl;
	CrackResidual2D crack_residual(30.f, 6, thread_number);
	crack_residual.prepare(poi_queue);
	crack_residual.compute(ref_img, tar_img, poi_queue);

	const Eigen::MatrixXf& residual = crack_residual.residualMap();

	auto meanResidualNear = [&](int x_center, int half_width) -> pair<float, int>
	{
		float sum = 0.f;
		int count = 0;
		for (int y = 40; y <= height - 40; y += 5)
		{
			for (int x = x_center - half_width; x <= x_center + half_width; x++)
			{
				float r = residual(y, x);
				if (r >= 0.f) { sum += r; count++; }
			}
		}
		return { count > 0 ? sum / count : -1.f, count };
	};

	auto near_line = meanResidualNear(x0, 3);
	auto far_left = meanResidualNear(x0 - 50, 3);
	auto far_right = meanResidualNear(x0 + 50, 3);

	cout << "  mean residual near the line (x0+/-3): " << near_line.first << " (" << near_line.second << " valid pixels)" << endl;
	cout << "  mean residual far left (x0-50+/-3):   " << far_left.first << " (" << far_left.second << " valid pixels)" << endl;
	cout << "  mean residual far right (x0+50+/-3):  " << far_right.first << " (" << far_right.second << " valid pixels)" << endl;

	bool residual_ok = near_line.second > 0 && far_left.second > 0 && far_right.second > 0
		&& near_line.first > 3.f * far_left.first
		&& near_line.first > 3.f * far_right.first;
	cout << "  " << (residual_ok ? "PASS" : "FAIL") << ": residual near the line is markedly higher than far from it, "
		"on a discontinuity every POI's own ZNCC missed" << endl;
	if (!residual_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
