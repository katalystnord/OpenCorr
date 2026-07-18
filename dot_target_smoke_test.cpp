/*
 Local build-verification test for CameraCalibrator's dot-target detection,
 fork issues #3/#10.

 Generates synthetic dot-target images from a KNOWN ground-truth camera
 model (multiple poses, perspective-warped from a flat dot-grid texture
 with 3 donut-shaped marker dots) and verifies:
   1. detectDots() (exercised indirectly through addImage/addImagePair)
      recovers the correct number of grid points and correctly locates the
      3 donut markers (checked by the calibration converging at all --
      wrong marker identification scrambles the origin/axis assignment and
      produces an obviously-wrong or non-converging solve).
   2. calibrate()/calibrateStereo() recover the ground-truth camera
      parameters within the same tolerance calibration_smoke_test.cpp uses
      for the checkerboard path.

 Same scope note as calibration_smoke_test.cpp: pure perspective warp, zero
 lens distortion in the ground truth.
*/

#include <iostream>
#include <random>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

//renders a synthetic dot-target image as seen by a pinhole camera: a
//num_fiducials_x x num_fiducials_y grid of filled black dots on white,
//dot_spacing units apart, with 3 dots replaced by donut markers (outer
//filled circle + inner hole) at the origin/x-axis-marker/y-axis-marker
//grid positions
static cv::Mat renderDotTarget(int num_fiducials_x, int num_fiducials_y, float dot_spacing,
	int origin_x, int origin_y, int origin_to_x_marker, int origin_to_y_marker,
	const cv::Mat& camera_matrix, const cv::Mat& rvec, const cv::Mat& tvec,
	int image_width, int image_height)
{
	const int px_per_dot = 300;
	const int margin_dots = 1;
	const int dot_radius = (int)(px_per_dot * 0.30f);
	const int hole_radius = (int)(px_per_dot * 0.13f);

	cv::Mat texture((num_fiducials_y + 2 * margin_dots) * px_per_dot, (num_fiducials_x + 2 * margin_dots) * px_per_dot,
		CV_8UC1, cv::Scalar(255));

	for (int r = 0; r < num_fiducials_y; r++)
	{
		for (int c = 0; c < num_fiducials_x; c++)
		{
			cv::Point center((c + margin_dots) * px_per_dot + px_per_dot / 2, (r + margin_dots) * px_per_dot + px_per_dot / 2);
			cv::circle(texture, center, dot_radius, cv::Scalar(0), cv::FILLED);

			bool is_origin = (c == origin_x && r == origin_y);
			bool is_x_marker = (c == origin_x + origin_to_x_marker - 1 && r == origin_y);
			bool is_y_marker = (c == origin_x && r == origin_y + origin_to_y_marker - 1);
			if (is_origin || is_x_marker || is_y_marker)
			{
				cv::circle(texture, center, hole_radius, cv::Scalar(255), cv::FILLED);
			}
		}
	}

	//map the ENTIRE texture onto the exact physical region it represents (same
	//full-bbox-mapping approach as calibration_smoke_test.cpp's renderCheckerboard():
	//a dot at grid (c,r) sits at texture pixel ((c+margin_dots)*px_per_dot + px_per_dot/2, ...)
	//and must map to world (c*dot_spacing, r*dot_spacing, 0)
	float edge = (margin_dots + 0.5f) * dot_spacing;
	vector<cv::Point3f> world_corners = {
		{-edge, -edge, 0.f},
		{(num_fiducials_x - 1) * dot_spacing + edge, -edge, 0.f},
		{(num_fiducials_x - 1) * dot_spacing + edge, (num_fiducials_y - 1) * dot_spacing + edge, 0.f},
		{-edge, (num_fiducials_y - 1) * dot_spacing + edge, 0.f}
	};
	vector<cv::Point2f> texture_corners = {
		{0.f, 0.f}, {(float)texture.cols, 0.f},
		{(float)texture.cols, (float)texture.rows}, {0.f, (float)texture.rows}
	};

	vector<cv::Point2f> projected;
	cv::Mat dist_zero = cv::Mat::zeros(1, 5, CV_64F);
	cv::projectPoints(world_corners, rvec, tvec, camera_matrix, dist_zero, projected);

	cv::Mat homography = cv::getPerspectiveTransform(texture_corners, projected);
	cv::Mat image;
	cv::warpPerspective(texture, image, homography, cv::Size(image_width, image_height),
		cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(180));

	return image;
}

int main()
{
	int failures = 0;

	int num_fiducials_x = 7, num_fiducials_y = 5;
	float dot_spacing = 25.f;
	int origin_x = 0, origin_y = 0;
	int origin_to_x_marker = num_fiducials_x; //x-axis marker at the far end of the row
	int origin_to_y_marker = num_fiducials_y; //y-axis marker at the far end of the column

	double gt_fx = 900.0, gt_fy = 900.0, gt_cx = 320.0, gt_cy = 240.0;
	cv::Mat gt_camera_matrix = (cv::Mat_<double>(3, 3) << gt_fx, 0, gt_cx, 0, gt_fy, gt_cy, 0, 0, 1);
	int image_width = 640, image_height = 480;

	std::mt19937 rng(7);
	std::uniform_real_distribution<float> tilt(-0.4f, 0.4f);
	std::uniform_real_distribution<float> in_plane_shift(-60.f, 60.f);

	//--- single-camera calibration ---
	cout << "=== Single-camera dot-target calibration ===" << endl;
	CameraCalibrator mono_calibrator(num_fiducials_x, num_fiducials_y, dot_spacing,
		origin_x, origin_y, origin_to_x_marker, origin_to_y_marker);

	int n_poses = 16;
	for (int i = 0; i < n_poses; i++)
	{
		float tx = in_plane_shift(rng), ty = in_plane_shift(rng), tz = 550.f + (i % 4) * 100.f;
		float rx = tilt(rng), ry = tilt(rng), rz = tilt(rng) * 0.5f;
		cv::Mat rvec = (cv::Mat_<double>(3, 1) << rx, ry, rz);
		cv::Mat tvec = (cv::Mat_<double>(3, 1) << tx - num_fiducials_x * dot_spacing / 2, ty - num_fiducials_y * dot_spacing / 2, tz);

		cv::Mat img = renderDotTarget(num_fiducials_x, num_fiducials_y, dot_spacing,
			origin_x, origin_y, origin_to_x_marker, origin_to_y_marker,
			gt_camera_matrix, rvec, tvec, image_width, image_height);
		string path = "/tmp/oc_dot_mono_" + to_string(i) + ".png";
		cv::imwrite(path, img);
		if (!mono_calibrator.addImage(path))
		{
			cout << "  warning: pose " << i << " not detected" << endl;
		}
	}

	cout << "  " << mono_calibrator.imageCount() << " / " << n_poses << " images used" << endl;

	Calibration mono_camera;
	bool mono_ok = false;
	float mono_rms = 0.f;
	try
	{
		mono_rms = mono_calibrator.calibrate(mono_camera);
		cout << "  RMS reprojection error: " << mono_rms << " px" << endl;
		cout << "  recovered fx=" << mono_camera.intrinsics.fx << " fy=" << mono_camera.intrinsics.fy
			<< " cx=" << mono_camera.intrinsics.cx << " cy=" << mono_camera.intrinsics.cy << endl;
		cout << "  ground truth fx=" << gt_fx << " fy=" << gt_fy << " cx=" << gt_cx << " cy=" << gt_cy << endl;

		mono_ok = mono_calibrator.imageCount() >= 8
			&& mono_rms < 1.5f
			&& fabs(mono_camera.intrinsics.fx - gt_fx) < 8.0
			&& fabs(mono_camera.intrinsics.fy - gt_fy) < 8.0
			&& fabs(mono_camera.intrinsics.cx - gt_cx) < 8.0
			&& fabs(mono_camera.intrinsics.cy - gt_cy) < 8.0;
	}
	catch (const std::string& e)
	{
		cout << "  exception: " << e << endl;
	}
	cout << "  " << (mono_ok ? "PASS" : "FAIL") << endl;
	if (!mono_ok) failures++;

	//--- stereo calibration ---
	cout << endl << "=== Stereo dot-target calibration ===" << endl;
	CameraCalibrator stereo_calibrator(num_fiducials_x, num_fiducials_y, dot_spacing,
		origin_x, origin_y, origin_to_x_marker, origin_to_y_marker);

	cv::Mat r_rel = cv::Mat::eye(3, 3, CV_64F);
	cv::Mat t_rel = (cv::Mat_<double>(3, 1) << 120.0, 0.0, 0.0);

	for (int i = 0; i < n_poses; i++)
	{
		float tx = in_plane_shift(rng), ty = in_plane_shift(rng), tz = 550.f + (i % 4) * 100.f;
		float rx = tilt(rng), ry = tilt(rng), rz = tilt(rng) * 0.5f;
		cv::Mat rvec_l = (cv::Mat_<double>(3, 1) << rx, ry, rz);
		cv::Mat tvec_l = (cv::Mat_<double>(3, 1) << tx - num_fiducials_x * dot_spacing / 2, ty - num_fiducials_y * dot_spacing / 2, tz);

		cv::Mat R_l;
		cv::Rodrigues(rvec_l, R_l);
		cv::Mat R_r = r_rel * R_l;
		cv::Mat T_r = r_rel * tvec_l + t_rel;
		cv::Mat rvec_r;
		cv::Rodrigues(R_r, rvec_r);

		cv::Mat img_l = renderDotTarget(num_fiducials_x, num_fiducials_y, dot_spacing,
			origin_x, origin_y, origin_to_x_marker, origin_to_y_marker,
			gt_camera_matrix, rvec_l, tvec_l, image_width, image_height);
		cv::Mat img_r = renderDotTarget(num_fiducials_x, num_fiducials_y, dot_spacing,
			origin_x, origin_y, origin_to_x_marker, origin_to_y_marker,
			gt_camera_matrix, rvec_r, T_r, image_width, image_height);

		string path_l = "/tmp/oc_dot_stereo_l_" + to_string(i) + ".png";
		string path_r = "/tmp/oc_dot_stereo_r_" + to_string(i) + ".png";
		cv::imwrite(path_l, img_l);
		cv::imwrite(path_r, img_r);

		if (!stereo_calibrator.addImagePair(path_l, path_r))
		{
			cout << "  warning: pose " << i << " not detected in both views" << endl;
		}
	}

	cout << "  " << stereo_calibrator.pairCount() << " / " << n_poses << " pairs used" << endl;

	bool stereo_ok = false;
	try
	{
		Calibration left_camera, right_camera;
		CameraExtrinsics right_extrinsics;
		float stereo_rms = stereo_calibrator.calibrateStereo(left_camera, right_camera, right_extrinsics);
		cout << "  RMS reprojection error: " << stereo_rms << " px" << endl;
		cout << "  recovered baseline: tx=" << right_extrinsics.tx << " ty=" << right_extrinsics.ty << " tz=" << right_extrinsics.tz
			<< " (ground truth tx=120, ty=0, tz=0)" << endl;

		stereo_ok = stereo_calibrator.pairCount() >= 8
			&& stereo_rms < 1.5f
			&& fabs(right_extrinsics.tx - 120.0) < 5.0
			&& fabs(right_extrinsics.ty) < 5.0
			&& fabs(right_extrinsics.tz) < 5.0;
	}
	catch (const std::string& e)
	{
		cout << "  exception: " << e << endl;
	}
	cout << "  " << (stereo_ok ? "PASS" : "FAIL") << endl;
	if (!stereo_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
