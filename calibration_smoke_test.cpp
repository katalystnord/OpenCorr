/*
 Local build-verification test for CameraCalibrator, issue #3.

 Generates synthetic checkerboard images from a KNOWN ground-truth camera
 model (multiple poses, perspective-warped from a flat checkerboard
 texture) and verifies the calibration solve recovers those parameters
 within tolerance, then verifies the stereo path's epipolar-residual
 metric behaves correctly on both a genuinely well-calibrated pair set and
 a deliberately-corrupted one (to confirm the metric actually detects bad
 calibration, not just that it produces *a* number).

 Scope note: images are generated via a pure perspective (homography) warp
 with zero lens distortion in the ground truth, so this exercises the
 calibration orchestration (corner detection -> calibrateCamera/
 stereoCalibrate -> parameter recovery -> epipolar quality metric) but not
 distortion-coefficient recovery specifically.
*/

#include <iostream>
#include <random>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

//renders a synthetic checkerboard image as seen by a pinhole camera with the
//given intrinsics, positioned/oriented by rvec/tvec relative to a flat
//checkerboard target lying in the z=0 plane (board_width x board_height
//internal corners, square_size units apart, with a 1-square quiet margin)
static cv::Mat renderCheckerboard(int board_width, int board_height, float square_size,
	const cv::Mat& camera_matrix, const cv::Mat& rvec, const cv::Mat& tvec,
	int image_width, int image_height)
{
	int squares_x = board_width + 1;
	int squares_y = board_height + 1;
	int px_per_square = 240; //high-res texture to keep the perspective-warped result anti-aliased at oblique angles

	//build the flat texture (with a 1-square quiet margin so all corners are interior)
	cv::Mat texture((squares_y + 2) * px_per_square, (squares_x + 2) * px_per_square, CV_8UC1, cv::Scalar(255));
	for (int r = 0; r < squares_y; r++)
		for (int c = 0; c < squares_x; c++)
			if ((r + c) % 2 == 0)
				cv::rectangle(texture,
					cv::Point((c + 1) * px_per_square, (r + 1) * px_per_square),
					cv::Point((c + 2) * px_per_square, (r + 2) * px_per_square),
					cv::Scalar(0), cv::FILLED);

	//map the ENTIRE texture (including both margins) onto the exact physical region it
	//represents -- picking a hand-selected sub-region here previously introduced a subtle
	//off-by-one-square mismatch between world_corners and texture_corners (225 vs 250
	//physical units, an ~11% stretch) that silently biased every recovered parameter by
	//the same ~11%, which is exactly what showed up in the baseline recovery below.
	vector<cv::Point3f> world_corners = {
		{-square_size, -square_size, 0.f}, {(squares_x + 1) * square_size, -square_size, 0.f},
		{(squares_x + 1) * square_size, (squares_y + 1) * square_size, 0.f}, {-square_size, (squares_y + 1) * square_size, 0.f}
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
		cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(128));

	return image;
}

int main()
{
	int failures = 0;

	int board_width = 8, board_height = 6;
	float square_size = 25.f; //arbitrary physical unit, e.g. mm

	//ground-truth intrinsics
	double gt_fx = 900.0, gt_fy = 900.0, gt_cx = 320.0, gt_cy = 240.0;
	cv::Mat gt_camera_matrix = (cv::Mat_<double>(3, 3) << gt_fx, 0, gt_cx, 0, gt_fy, gt_cy, 0, 0, 1);
	int image_width = 640, image_height = 480;

	//--- single-camera calibration ---
	cout << "=== Single-camera calibration ===" << endl;
	CameraCalibrator mono_calibrator(board_width, board_height, square_size);

	//pose diversity matters a lot for conditioning: a set of near-frontal, near-identical-depth
	//views leaves focal length weakly separable from depth (a known calibration ill-conditioning
	//issue, not specific to this implementation). Real calibration guides recommend tilting the
	//target 15-45 degrees across shots for exactly this reason (matches the VIC-Snap/GOM/Istra4D
	//guidance from the DIC competitive research) -- so this test does the same.
	std::mt19937 rng(42);
	std::uniform_real_distribution<float> tilt(-0.5f, 0.5f); //~ +/-29 degrees
	std::uniform_real_distribution<float> in_plane_shift(-80.f, 80.f);

	int n_poses = 16;
	for (int i = 0; i < n_poses; i++)
	{
		float tx = in_plane_shift(rng), ty = in_plane_shift(rng), tz = 500.f + (i % 4) * 120.f;
		float rx = tilt(rng), ry = tilt(rng), rz = tilt(rng) * 0.5f;
		cv::Mat rvec = (cv::Mat_<double>(3, 1) << rx, ry, rz);
		cv::Mat tvec = (cv::Mat_<double>(3, 1) << tx - board_width * square_size / 2, ty - board_height * square_size / 2, tz);

		cv::Mat img = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec, tvec, image_width, image_height);
		string path = "/tmp/oc_cal_mono_" + to_string(i) + ".png";
		cv::imwrite(path, img);
		if (!mono_calibrator.addImage(path))
		{
			cout << "  warning: pose " << i << " not detected" << endl;
		}
	}

	cout << "  " << mono_calibrator.imageCount() << " / " << n_poses << " images used" << endl;

	Calibration mono_camera;
	float mono_rms = mono_calibrator.calibrate(mono_camera);
	cout << "  RMS reprojection error: " << mono_rms << " px" << endl;
	cout << "  recovered fx=" << mono_camera.intrinsics.fx << " fy=" << mono_camera.intrinsics.fy
		<< " cx=" << mono_camera.intrinsics.cx << " cy=" << mono_camera.intrinsics.cy << endl;
	cout << "  ground truth fx=" << gt_fx << " fy=" << gt_fy << " cx=" << gt_cx << " cy=" << gt_cy << endl;

	//calibrate() must leave the Calibration ready for Stereovision/EpipolarSearch to use
	//directly (i.e. it must call prepare() itself) -- every other caller of Calibration in
	//this codebase calls prepare() before use, and map_x/map_y (populated by prepare(),
	//oc_calibration.h) staying empty is exactly the bug this fix corrects
	bool mono_prepared = mono_camera.map_x.rows() == image_height && mono_camera.map_x.cols() == image_width;

	bool mono_ok = mono_calibrator.imageCount() >= 8
		&& mono_rms < 1.0f
		&& fabs(mono_camera.intrinsics.fx - gt_fx) < 5.0
		&& fabs(mono_camera.intrinsics.fy - gt_fy) < 5.0
		&& fabs(mono_camera.intrinsics.cx - gt_cx) < 5.0
		&& fabs(mono_camera.intrinsics.cy - gt_cy) < 5.0
		&& mono_prepared;
	cout << "  calibrate() left the Calibration prepare()'d (map_x is " << mono_camera.map_x.rows()
		<< "x" << mono_camera.map_x.cols() << ", expected " << image_height << "x" << image_width << "): "
		<< (mono_prepared ? "yes" : "no") << endl;
	cout << "  " << (mono_ok ? "PASS" : "FAIL") << endl;
	if (!mono_ok) failures++;

	//--- stereo calibration ---
	cout << endl << "=== Stereo calibration ===" << endl;
	CameraCalibrator stereo_calibrator(board_width, board_height, square_size);

	//right camera offset 120 units along x from the left camera, no relative rotation
	cv::Mat r_rel = cv::Mat::eye(3, 3, CV_64F);
	cv::Mat t_rel = (cv::Mat_<double>(3, 1) << 120.0, 0.0, 0.0);

	for (int i = 0; i < n_poses; i++)
	{
		float tx = in_plane_shift(rng), ty = in_plane_shift(rng), tz = 500.f + (i % 4) * 120.f;
		float rx = tilt(rng), ry = tilt(rng), rz = tilt(rng) * 0.5f;
		cv::Mat rvec_l = (cv::Mat_<double>(3, 1) << rx, ry, rz);
		cv::Mat tvec_l = (cv::Mat_<double>(3, 1) << tx - board_width * square_size / 2, ty - board_height * square_size / 2, tz);

		//right camera sees the same board, from a pose offset by the fixed relative transform
		cv::Mat R_l;
		cv::Rodrigues(rvec_l, R_l);
		cv::Mat R_r = r_rel * R_l;
		cv::Mat T_r = r_rel * tvec_l + t_rel;
		cv::Mat rvec_r;
		cv::Rodrigues(R_r, rvec_r);

		cv::Mat img_l = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec_l, tvec_l, image_width, image_height);
		cv::Mat img_r = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec_r, T_r, image_width, image_height);

		string path_l = "/tmp/oc_cal_stereo_l_" + to_string(i) + ".png";
		string path_r = "/tmp/oc_cal_stereo_r_" + to_string(i) + ".png";
		cv::imwrite(path_l, img_l);
		cv::imwrite(path_r, img_r);

		if (!stereo_calibrator.addImagePair(path_l, path_r))
		{
			cout << "  warning: pose " << i << " not detected in both views" << endl;
		}
	}

	cout << "  " << stereo_calibrator.pairCount() << " / " << n_poses << " pairs used" << endl;

	Calibration left_camera, right_camera;
	CameraExtrinsics right_extrinsics;
	float stereo_rms = stereo_calibrator.calibrateStereo(left_camera, right_camera, right_extrinsics);
	cout << "  RMS reprojection error: " << stereo_rms << " px" << endl;
	cout << "  recovered baseline: tx=" << right_extrinsics.tx << " ty=" << right_extrinsics.ty << " tz=" << right_extrinsics.tz
		<< " (ground truth tx=120, ty=0, tz=0)" << endl;

	float mean_epi = 0.f;
	for (float e : stereo_calibrator.epipolarResiduals()) mean_epi += e;
	mean_epi /= stereo_calibrator.epipolarResiduals().size();
	cout << "  mean epipolar residual: " << mean_epi << endl;

	bool stereo_prepared = left_camera.map_x.rows() == image_height && left_camera.map_x.cols() == image_width
		&& right_camera.map_x.rows() == image_height && right_camera.map_x.cols() == image_width;
	cout << "  calibrateStereo() left both Calibrations prepare()'d: " << (stereo_prepared ? "yes" : "no") << endl;

	bool stereo_ok = stereo_calibrator.pairCount() >= 8
		&& stereo_rms < 1.0f
		&& fabs(right_extrinsics.tx - 120.0) < 3.0
		&& fabs(right_extrinsics.ty) < 3.0
		&& fabs(right_extrinsics.tz) < 3.0
		&& mean_epi < 0.5f
		&& stereo_prepared;
	cout << "  " << (stereo_ok ? "PASS" : "FAIL") << endl;
	if (!stereo_ok) failures++;

	//--- sanity check: does the epipolar residual actually detect a BAD calibration? ---
	cout << endl << "=== Epipolar residual sensitivity check ===" << endl;
	cout << "  (deliberately mismatched pairs -- right image paired with the WRONG left image,"
		<< endl << "   simulating a miscalibrated/misidentified stereo rig)" << endl;
	CameraCalibrator bad_calibrator(board_width, board_height, square_size);
	for (int i = 0; i < n_poses; i++)
	{
		string path_l = "/tmp/oc_cal_stereo_l_" + to_string(i) + ".png";
		string path_r = "/tmp/oc_cal_stereo_r_" + to_string((i + 1) % n_poses) + ".png"; //off-by-one mismatch
		bad_calibrator.addImagePair(path_l, path_r);
	}
	Calibration bad_left, bad_right;
	CameraExtrinsics bad_extrinsics;
	bad_calibrator.calibrateStereo(bad_left, bad_right, bad_extrinsics);
	float bad_mean_epi = 0.f;
	for (float e : bad_calibrator.epipolarResiduals()) bad_mean_epi += e;
	bad_mean_epi /= bad_calibrator.epipolarResiduals().size();
	cout << "  mean epipolar residual (mismatched pairs): " << bad_mean_epi
		<< " vs. correctly-paired: " << mean_epi << endl;

	bool sensitivity_ok = bad_mean_epi > mean_epi * 5.f; //should be dramatically worse, not just noisier
	cout << "  " << (sensitivity_ok ? "PASS" : "FAIL") << " (metric correctly distinguishes good from bad calibration)" << endl;
	if (!sensitivity_ok) failures++;

	//--- epipolar residual must be evaluated in the UNDISTORTED domain ---
	//this doesn't exercise CameraCalibrator directly (the rest of this file's synthetic
	//images use zero ground-truth distortion, by design -- see the file's own scope note),
	//it isolates the exact principle epipolarResiduals() relies on: F only satisfies the
	//epipolar constraint for undistorted points, so feeding it raw distorted points (the bug
	//this fix corrects) should inflate the residual well above the undistorted case.
	cout << endl << "=== Epipolar residual must use the undistorted domain, not raw distorted points ===" << endl;
	cv::Mat cam_l = (cv::Mat_<double>(3, 3) << 900, 0, 320, 0, 900, 240, 0, 0, 1);
	cv::Mat cam_r = cam_l.clone();
	cv::Mat dist = (cv::Mat_<double>(1, 5) << -0.35, 0.15, 0, 0, 0); //real, non-negligible radial distortion
	cv::Mat rvec_stub = cv::Mat::zeros(3, 1, CV_64F);
	cv::Mat tvec_l_stub = cv::Mat::zeros(3, 1, CV_64F);
	cv::Mat tvec_r_stub = (cv::Mat_<double>(3, 1) << 120.0, 0.0, 0.0);

	vector<cv::Point3f> world_pts;
	for (float x = -30.f; x <= 30.f; x += 10.f)
		for (float y = -30.f; y <= 30.f; y += 10.f)
			world_pts.push_back(cv::Point3f(x, y, 400.f)); //well in front of both cameras

	vector<cv::Point2f> ideal_l, ideal_r, distorted_l, distorted_r;
	cv::Mat dist_zero5 = cv::Mat::zeros(1, 5, CV_64F);
	cv::projectPoints(world_pts, rvec_stub, tvec_l_stub, cam_l, dist_zero5, ideal_l);
	cv::projectPoints(world_pts, rvec_stub, tvec_r_stub, cam_r, dist_zero5, ideal_r);
	cv::projectPoints(world_pts, rvec_stub, tvec_l_stub, cam_l, dist, distorted_l);
	cv::projectPoints(world_pts, rvec_stub, tvec_r_stub, cam_r, dist, distorted_r);

	//F fit from the IDEAL (undistorted) correspondences -- exactly what cv::stereoCalibrate
	//itself returns, since it always operates in the undistorted/normalized domain internally
	cv::Mat F_synth = cv::findFundamentalMat(ideal_l, ideal_r, cv::FM_8POINT);

	auto meanEpipolarResidual = [&](const vector<cv::Point2f>& pl, const vector<cv::Point2f>& pr) -> double
	{
		vector<cv::Vec3f> lines_l, lines_r;
		cv::computeCorrespondEpilines(pl, 1, F_synth, lines_l);
		cv::computeCorrespondEpilines(pr, 2, F_synth, lines_r);
		double total = 0.0;
		for (size_t j = 0; j < pl.size(); j++)
		{
			total += std::fabs(pl[j].x * lines_r[j][0] + pl[j].y * lines_r[j][1] + lines_r[j][2])
				+ std::fabs(pr[j].x * lines_l[j][0] + pr[j].y * lines_l[j][1] + lines_l[j][2]);
		}
		return total / pl.size();
	};

	double residual_undistorted = meanEpipolarResidual(ideal_l, ideal_r);
	double residual_distorted = meanEpipolarResidual(distorted_l, distorted_r);
	cout << "  mean residual, undistorted points vs F: " << residual_undistorted << endl;
	cout << "  mean residual, RAW DISTORTED points vs the same F: " << residual_distorted << endl;

	bool undistort_domain_ok = residual_undistorted < 0.05 && residual_distorted > residual_undistorted * 20.0;
	cout << "  " << (undistort_domain_ok ? "PASS" : "FAIL")
		<< ": undistorted points satisfy F almost exactly, raw distorted points don't -- "
		<< "confirms epipolarResiduals() must undistort before evaluating, matching the actual fix" << endl;
	if (!undistort_domain_ok) failures++;

	//--- regression: a genuinely mismatched-resolution stereo rig must not have its
	//right-camera images spuriously rejected as a "size mismatch" against the left
	//camera's own resolution (image_size used to be a single value shared across both) ---
	cout << endl << "=== Mismatched left/right resolutions must both be accepted ===" << endl;
	{
		CameraCalibrator mismatched_calibrator(board_width, board_height, square_size);
		int right_width = 800, right_height = 600; //deliberately different from the left's 640x480

		cv::Mat rvec_l = (cv::Mat_<double>(3, 1) << 0.1, 0.05, 0.0);
		cv::Mat tvec_l = (cv::Mat_<double>(3, 1) << -board_width * square_size / 2, -board_height * square_size / 2, 600.0);
		cv::Mat R_l;
		cv::Rodrigues(rvec_l, R_l);
		cv::Mat R_r = r_rel * R_l;
		cv::Mat T_r = r_rel * tvec_l + t_rel;
		cv::Mat rvec_r;
		cv::Rodrigues(R_r, rvec_r);

		//left camera's own intrinsics reused for the right image too -- this test is only
		//about whether detection/size-tracking accepts the pair, not stereo parameter
		//recovery, so a physically-accurate right-camera intrinsic matrix isn't needed
		cv::Mat img_l = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec_l, tvec_l, image_width, image_height);
		cv::Mat img_r = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec_r, T_r, right_width, right_height);

		string path_l = "/tmp/oc_cal_mismatched_l.png";
		string path_r = "/tmp/oc_cal_mismatched_r.png";
		cv::imwrite(path_l, img_l);
		cv::imwrite(path_r, img_r);

		bool pair_added = mismatched_calibrator.addImagePair(path_l, path_r);
		cout << "  left " << image_width << "x" << image_height << ", right " << right_width << "x" << right_height
			<< " -- pair accepted: " << pair_added << endl;
		cout << "  " << (pair_added ? "PASS" : "FAIL")
			<< ": right camera's own (different) resolution isn't rejected against the left's" << endl;
		if (!pair_added) failures++;
	}

	//--- regression: calibrating from fewer than 2 views must throw a clear error, not
	//silently hand cv::calibrateCamera/cv::stereoCalibrate a mathematically underdetermined
	//problem that can numerically "succeed" while producing a degenerate/near-singular result ---
	cout << endl << "=== Fewer than 2 views must throw, not silently under-solve ===" << endl;
	{
		CameraCalibrator one_view_calibrator(board_width, board_height, square_size);
		cv::Mat rvec = (cv::Mat_<double>(3, 1) << 0.1, 0.05, 0.0);
		cv::Mat tvec = (cv::Mat_<double>(3, 1) << -board_width * square_size / 2, -board_height * square_size / 2, 600.0);
		cv::Mat img = renderCheckerboard(board_width, board_height, square_size, gt_camera_matrix, rvec, tvec, image_width, image_height);
		string path = "/tmp/oc_cal_one_view.png";
		cv::imwrite(path, img);
		one_view_calibrator.addImage(path);

		bool threw = false;
		try
		{
			Calibration dummy_camera;
			one_view_calibrator.calibrate(dummy_camera);
		}
		catch (std::string&)
		{
			threw = true;
		}
		cout << "  1 image, calibrate() threw: " << threw << endl;
		cout << "  " << (threw ? "PASS" : "FAIL") << ": calibrate() rejects a single view instead of under-solving silently" << endl;
		if (!threw) failures++;

		CameraCalibrator one_pair_calibrator(board_width, board_height, square_size);
		one_pair_calibrator.addImagePair("/tmp/oc_cal_mismatched_l.png", "/tmp/oc_cal_mismatched_r.png");
		bool stereo_threw = false;
		try
		{
			Calibration dummy_left, dummy_right;
			CameraExtrinsics dummy_extrinsics;
			one_pair_calibrator.calibrateStereo(dummy_left, dummy_right, dummy_extrinsics);
		}
		catch (std::string&)
		{
			stereo_threw = true;
		}
		cout << "  1 pair, calibrateStereo() threw: " << stereo_threw << endl;
		cout << "  " << (stereo_threw ? "PASS" : "FAIL") << ": calibrateStereo() rejects a single pair instead of under-solving silently" << endl;
		if (!stereo_threw) failures++;
	}

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
