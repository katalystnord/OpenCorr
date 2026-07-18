/*
 Local build-verification test for PhaseCorrelation2D, issue #6.

 1. Ground truth: circularly shift a real image by a KNOWN integer pixel
    offset (a circular/toroidal shift is what FFT-based phase correlation
    natively assumes -- no edge artifacts to confound the result) and
    verify the recovered displacement is exact.
 2. Real motion: cross-check the coarse whole-image estimate against
    FFTCC2D's own subset-level result on the same well-textured image pair
    already used elsewhere in this repo's examples.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;
	omp_set_num_threads(4); //must match the thread_number FFTCC2D below is constructed with

	//--- check 1: exact recovery of a known circular shift ---
	cout << "=== Ground truth: circularly-shifted image ===" << endl;
	Image2D ref_img("examples/2d_dic/oht_cfrp_0.bmp");
	int w = ref_img.width, h = ref_img.height;

	int shift_x = 17, shift_y = -9;
	Image2D shifted_img(w, h);
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			int src_x = ((x - shift_x) % w + w) % w;
			int src_y = ((y - shift_y) % h + h) % h;
			shifted_img.cv_mat.at<uchar>(y, x) = ref_img.cv_mat.at<uchar>(src_y, src_x);
		}
	}
	cv::cv2eigen(shifted_img.cv_mat, shifted_img.eg_mat);

	PhaseCorrelation2D phase_corr(w, h);
	float u = 0.f, v = 0.f;
	float confidence = phase_corr.compute(ref_img, shifted_img, u, v);

	cout << "  injected shift: (" << shift_x << ", " << shift_y << ")" << endl;
	cout << "  recovered:      (" << u << ", " << v << "), peak magnitude=" << confidence << endl;
	bool exact_ok = (int)u == shift_x && (int)v == shift_y;
	cout << "  " << (exact_ok ? "PASS" : "FAIL") << ": exact recovery of circular shift" << endl;
	if (!exact_ok) failures++;

	//--- check 2: cross-check against FFTCC2D on real motion ---
	cout << endl << "=== Cross-check against FFTCC2D on real motion ===" << endl;
	Image2D tar_img("examples/2d_dic/oht_cfrp_4.bmp");

	PhaseCorrelation2D phase_corr2(w, h);
	float u2 = 0.f, v2 = 0.f;
	float confidence2 = phase_corr2.compute(ref_img, tar_img, u2, v2);
	cout << "  whole-image phase-correlation estimate: (" << u2 << ", " << v2 << "), peak magnitude=" << confidence2 << endl;

	//sample FFTCC2D at a grid of POIs and report the median subset-level displacement,
	//for comparison -- the two shouldn't need to match exactly (phase correlation is a
	//single global rigid estimate; FFTCC2D is per-subset and can pick up local
	//deformation), but they should be in the same ballpark for this real dataset, which
	//is a small-deformation open-hole tension specimen, not a large rigid motion
	vector<POI2D> poi_queue;
	for (int y = 30; y < h - 30; y += 20)
		for (int x = 30; x < w - 30; x += 40)
			poi_queue.push_back(POI2D(Point2D((float)x, (float)y)));

	FFTCC2D fftcc(16, 16, 4);
	fftcc.setImages(ref_img, tar_img);
	fftcc.compute(poi_queue);

	vector<float> us, vs;
	for (auto& poi : poi_queue) { us.push_back(poi.deformation.u); vs.push_back(poi.deformation.v); }
	sort(us.begin(), us.end());
	sort(vs.begin(), vs.end());
	float median_u = us[us.size() / 2], median_v = vs[vs.size() / 2];
	cout << "  FFTCC2D median subset-level displacement: (" << median_u << ", " << median_v << ")" << endl;

	float diff = sqrt(pow(u2 - median_u, 2) + pow(v2 - median_v, 2));
	cout << "  difference: " << diff << " px" << endl;
	//phase correlation only resolves to integer-pixel precision by construction (it finds
	//a peak FFT bin, no sub-pixel refinement is implemented), so agreement within ~1.5px
	//of the subset-level median is the right bar here, not sub-pixel agreement
	bool crosscheck_ok = diff < 1.5f;
	cout << "  " << (crosscheck_ok ? "PASS" : "FAIL") << ": whole-image and subset-level estimates agree within integer-pixel tolerance" << endl;
	if (!crosscheck_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
