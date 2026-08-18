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

//builds an image with NO speckled region at all: a uniform field plus faint
//sensor noise. There is nothing here to segment, and the correct output is
//"no region", not a plausible-looking one.
static Image2D renderBlankFrame(int image_width, int image_height, double sigma, unsigned seed)
{
	Image2D image(image_width, image_height);

	std::mt19937 rng(seed);
	std::normal_distribution<double> noise(0.0, sigma);
	for (int r = 0; r < image_height; r++)
	{
		for (int c = 0; c < image_width; c++)
		{
			int v = (int)std::lround(128.0 + noise(rng));
			image.cv_mat.at<uchar>(r, c) = (uchar)std::max(0, std::min(255, v));
		}
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
	//detect() returns the Shape2D base (issue #17) since a future hole-aware detection
	//could return a different concrete type -- today it always actually constructs a
	//Polygon2D, so dynamic_cast back to it for numVertices()/rectangleIoU() below is safe
	//and also doubles as a check that this test still matches what detect() actually returns
	std::unique_ptr<Shape2D> detected = auto_roi.detect(image);
	Polygon2D* detected_polygon = dynamic_cast<Polygon2D*>(detected.get());

	bool roi_ok = false;
	if (detected_polygon != nullptr)
	{
		float iou = rectangleIoU(*detected_polygon, rect_x, rect_y, rect_w, rect_h);
		cout << "  detected polygon: " << detected_polygon->numVertices() << " vertices, bbox ["
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
	//This image is a coupon speckled across essentially the whole frame, so it
	//has no speckle-versus-background structure to segment. It is kept as a
	//check that real data does not crash the detector, and that whatever it
	//decides is HONEST -- either a region that is genuinely smaller than the
	//frame, or no region at all.
	//
	//The expectation here used to be "finds a plausible non-degenerate region",
	//which the old code satisfied by returning a bounding box covering 99.5% of
	//the image. A region identical to the frame is not a segmentation; it is a
	//refusal wearing a polygon, and asserting it as a success is how the
	//separability guard's failure stayed invisible.
	Image2D real_image("examples/2d_dic/oht_cfrp_0.bmp");
	AutoROI real_roi(20);
	std::unique_ptr<Shape2D> real_detected_polygon = real_roi.detect(real_image);

	bool real_ok = true;
	if (real_detected_polygon != nullptr)
	{
		int w = real_detected_polygon->getMaxX() - real_detected_polygon->getMinX() + 1;
		int h = real_detected_polygon->getMaxY() - real_detected_polygon->getMinY() + 1;
		double area_fraction = (double)w * h / ((double)real_image.width * real_image.height);
		cout << "  region " << w << " x " << h
			<< ", bbox area / image area = " << area_fraction << endl;
		//A "region" that is the whole frame tells the caller nothing.
		real_ok = area_fraction < 0.95;
		if (!real_ok) cout << "    that is the whole frame, which is not a segmentation" << endl;
	}
	else
	{
		cout << "  declined -- no speckle-versus-background structure in this frame" << endl;
	}
	cout << "  " << (real_ok ? "PASS" : "FAIL")
		<< ": real data yields either a genuine sub-region or an honest refusal" << endl;
	if (!real_ok) failures++;

	//--- a frame with nothing in it must be DECLINED, not segmented ---
	//
	//This is the case detect() documents itself as guarding against, and the
	//guard did not work. cv::normalize stretches pure sensor noise to full
	//contrast before thresholding, so a blank frame arrives at Otsu looking
	//exactly like a high-contrast one; the separability floor was supposed to
	//catch that and cannot, because Otsu separability answers "how well can
	//this be cut in two", which is ~0.64 for ANY unimodal spread -- four times
	//the 0.15 floor -- and independent of the noise amplitude, since
	//normalisation removes the scale.
	//
	//Found from the consumer side: a SurView test asserting that auto-detect
	//declines on a blank frame failed, and returned a 25-corner region covering
	//about 80% of the image.
	cout << endl << "=== A blank frame must be declined, not segmented ===" << endl;
	int blank_regions = 0;
	for (double sigma : { 0.5, 1.2, 3.0, 10.0 })
	{
		Image2D blank = renderBlankFrame(240, 200, sigma, 7);
		AutoROI blank_roi(15);
		std::unique_ptr<Shape2D> blank_detected = blank_roi.detect(blank);
		bool declined = (blank_detected == nullptr);
		if (!declined)
		{
			blank_regions++;
			cout << "    sigma " << sigma << ": returned a region of "
				<< (blank_detected->getMaxX() - blank_detected->getMinX() + 1) << " x "
				<< (blank_detected->getMaxY() - blank_detected->getMinY() + 1) << " px" << endl;
		}
		else
		{
			cout << "    sigma " << sigma << ": declined" << endl;
		}
	}
	cout << "  " << (blank_regions == 0 ? "PASS" : "FAIL")
		<< ": noise-only frames yield no region (" << blank_regions << " of 4 wrongly segmented)" << endl;
	if (blank_regions != 0) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
