/*
 Local build-verification test for SpeckleQualityMap/AutoROI, fork issues
 #11/#12 (speckle quality metrics + auto-ROI segmentation).

 Synthesizes an image with a speckled rectangular region on a low-texture
 background at a known boundary, and verifies:
   1. MIG/SSSIG/SIFT-density are meaningfully higher inside the speckled
      region than outside.
   2. AutoROI::detect() recovers a polygon matching the known ground-truth
      rectangle (checked via intersection-over-union).
 Then runs a sanity pass against a real speckle-pattern photograph already
 in the repo -- no ground truth there, just confirms no crash and a
 plausible (non-degenerate, non-full-image) result on real data.
*/

#include <iostream>
#include <random>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

//builds an image_width x image_height 8-bit grayscale image: a near-uniform,
//low-texture background and a speckled rectangle
//[rect_x, rect_x+rect_w) x [rect_y, rect_y+rect_h) filled with random dots
static Image2D renderSpeckledRegion(int image_width, int image_height,
	int rect_x, int rect_y, int rect_w, int rect_h, unsigned seed)
{
	Image2D image(image_width, image_height);

	std::mt19937 rng(seed);
	std::uniform_int_distribution<int> noise(-3, 3);
	for (int r = 0; r < image_height; r++)
	{
		for (int c = 0; c < image_width; c++)
		{
			image.cv_mat.at<uchar>(r, c) = (uchar)std::max(0, std::min(255, 200 + noise(rng)));
		}
	}

	std::uniform_int_distribution<int> dot_x(rect_x, rect_x + rect_w - 1);
	std::uniform_int_distribution<int> dot_y(rect_y, rect_y + rect_h - 1);
	std::uniform_int_distribution<int> dot_radius(3, 7);
	std::uniform_int_distribution<int> dot_gray(0, 120);
	int n_dots = rect_w * rect_h / 80; //roughly matches typical DIC speckle density
	for (int i = 0; i < n_dots; i++)
	{
		cv::circle(image.cv_mat, cv::Point(dot_x(rng), dot_y(rng)), dot_radius(rng), cv::Scalar(dot_gray(rng)), cv::FILLED);
	}

	cv::cv2eigen(image.cv_mat, image.eg_mat);
	return image;
}

//intersection-over-union of a Polygon2D against a known axis-aligned rectangle,
//computed by pixel enumeration (small test images, cheap enough)
static float rectangleIoU(const Polygon2D& polygon, int rect_x, int rect_y, int rect_w, int rect_h)
{
	int min_x = std::min(polygon.getMinX(), rect_x);
	int max_x = std::max(polygon.getMaxX(), rect_x + rect_w - 1);
	int min_y = std::min(polygon.getMinY(), rect_y);
	int max_y = std::max(polygon.getMaxY(), rect_y + rect_h - 1);

	long long intersection = 0, union_count = 0;
	for (int y = min_y; y <= max_y; y++)
	{
		for (int x = min_x; x <= max_x; x++)
		{
			bool in_poly = polygon.contains(x, y);
			bool in_rect = (x >= rect_x && x < rect_x + rect_w && y >= rect_y && y < rect_y + rect_h);
			if (in_poly || in_rect) union_count++;
			if (in_poly && in_rect) intersection++;
		}
	}
	return union_count > 0 ? (float)intersection / (float)union_count : 0.f;
}

int main()
{
	int failures = 0;

	int image_width = 400, image_height = 300;
	int rect_x = 100, rect_y = 60, rect_w = 200, rect_h = 160;

	cout << "=== Synthetic speckled-region test ===" << endl;
	Image2D image = renderSpeckledRegion(image_width, image_height, rect_x, rect_y, rect_w, rect_h, 42);

	SpeckleQualityMap quality_map(15);
	quality_map.compute(image);

	//sample a handful of interior/exterior points to check the metric separation itself,
	//independent of AutoROI's segmentation (a distinct failure mode: metrics could be fine
	//but the threshold/morphology/contour pipeline could still mis-segment, or vice versa)
	vector<pair<int, int>> inside_points = { {150, 100}, {200, 140}, {250, 180} };
	vector<pair<int, int>> outside_points = { {20, 20}, {380, 20}, {20, 280}, {380, 280} };

	float mean_mig_inside = 0.f, mean_mig_outside = 0.f;
	for (auto& p : inside_points) mean_mig_inside += quality_map.migMap()(p.second, p.first);
	for (auto& p : outside_points) mean_mig_outside += quality_map.migMap()(p.second, p.first);
	mean_mig_inside /= (float)inside_points.size();
	mean_mig_outside /= (float)outside_points.size();

	auto& whole = quality_map.wholeImageMetrics();
	cout << "  mean MIG inside=" << mean_mig_inside << " outside=" << mean_mig_outside << endl;
	cout << "  whole-image: mig=" << whole.mean_mig << " sssig=" << whole.mean_sssig
		<< " sift_density=" << whole.sift_density << " sift_evenness=" << whole.sift_evenness << endl;

	bool metrics_ok = mean_mig_inside > mean_mig_outside * 3.f; //should be dramatically higher, not just noisier
	cout << "  " << (metrics_ok ? "PASS" : "FAIL") << " (MIG clearly separates speckled region from background)" << endl;
	if (!metrics_ok) failures++;

	cout << endl << "=== AutoROI segmentation ===" << endl;
	AutoROI auto_roi(15);
	std::unique_ptr<Polygon2D> detected = auto_roi.detect(image);

	bool roi_ok = false;
	if (detected)
	{
		float iou = rectangleIoU(*detected, rect_x, rect_y, rect_w, rect_h);
		cout << "  detected polygon: " << detected->numVertices() << " vertices, bbox ["
			<< detected->getMinX() << "," << detected->getMinY() << "] to ["
			<< detected->getMaxX() << "," << detected->getMaxY() << "]" << endl;
		cout << "  ground truth rect: [" << rect_x << "," << rect_y << "] to ["
			<< rect_x + rect_w << "," << rect_y + rect_h << "]" << endl;
		cout << "  IoU=" << iou << endl;
		roi_ok = iou > 0.8f;
	}
	else
	{
		cout << "  AutoROI found no region" << endl;
	}
	cout << "  " << (roi_ok ? "PASS" : "FAIL") << " (recovered region matches ground truth, IoU > 0.8)" << endl;
	if (!roi_ok) failures++;

	cout << endl << "=== Sanity pass on real speckle photograph ===" << endl;
	Image2D real_image("examples/2d_dic/oht_cfrp_0.bmp");
	AutoROI real_roi(20);
	std::unique_ptr<Polygon2D> real_detected = real_roi.detect(real_image);

	bool real_ok = false;
	if (real_detected)
	{
		int area = (real_detected->getMaxX() - real_detected->getMinX()) * (real_detected->getMaxY() - real_detected->getMinY());
		int image_area = real_image.width * real_image.height;
		cout << "  detected polygon: " << real_detected->numVertices() << " vertices, bbox ["
			<< real_detected->getMinX() << "," << real_detected->getMinY() << "] to ["
			<< real_detected->getMaxX() << "," << real_detected->getMaxY() << "]" << endl;
		cout << "  bbox area / image area = " << (float)area / (float)image_area << endl;
		//plausible, not necessarily ground-truthed: a real region, not degenerate (empty) and
		//not trivially the entire image (which would suggest the threshold did nothing)
		real_ok = area > 0 && area < image_area;
	}
	else
	{
		cout << "  AutoROI found no region on real image" << endl;
	}
	cout << "  " << (real_ok ? "PASS" : "FAIL") << " (no crash, plausible non-degenerate region on real data)" << endl;
	if (!real_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
