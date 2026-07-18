/*
 * This file is part of OpenCorr, an open source C++ library for
 * study and development of 2D, 3D/stereo and volumetric
 * digital image correlation.
 *
 * Copyright (C) 2021-2025, Zhenyu Jiang <zhenyujiang@scut.edu.cn>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one from http://mozilla.org/MPL/2.0/.
 *
 * More information about OpenCorr can be found at https://www.opencorr.org/
 */

#include "oc_camera_calibrator.h"

namespace opencorr
{
	CameraCalibrator::CameraCalibrator(int board_width, int board_height, float square_size)
		: board_width(board_width), board_height(board_height), square_size(square_size), image_size_set(false)
	{
		for (int r = 0; r < board_height; r++)
		{
			for (int c = 0; c < board_width; c++)
			{
				object_points_template.push_back(cv::Point3f(c * square_size, r * square_size, 0.f));
			}
		}
	}

	CameraCalibrator::~CameraCalibrator() {}

	bool CameraCalibrator::detectCorners(const std::string& image_path, std::vector<cv::Point2f>& corners)
	{
		cv::Mat img = cv::imread(image_path, cv::IMREAD_GRAYSCALE);
		if (img.empty())
		{
			std::cerr << "CameraCalibrator: failed to load image: " << image_path << std::endl;
			return false;
		}

		if (!image_size_set)
		{
			image_size = img.size();
			image_size_set = true;
		}
		else if (img.size() != image_size)
		{
			std::cerr << "CameraCalibrator: image size mismatch, excluding: " << image_path << std::endl;
			return false;
		}

		bool found = cv::findChessboardCorners(img, cv::Size(board_width, board_height), corners,
			cv::CALIB_CB_ADAPTIVE_THRESH | cv::CALIB_CB_NORMALIZE_IMAGE);

		if (!found)
		{
			std::cerr << "CameraCalibrator: checkerboard not found, excluding: " << image_path << std::endl;
			return false;
		}

		cv::cornerSubPix(img, corners, cv::Size(11, 11), cv::Size(-1, -1),
			cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 30, 0.01));

		return true;
	}

	bool CameraCalibrator::addImage(const std::string& image_path)
	{
		std::vector<cv::Point2f> corners;
		if (!detectCorners(image_path, corners))
		{
			return false;
		}

		mono_image_points.push_back(corners);
		return true;
	}

	int CameraCalibrator::imageCount() const
	{
		return (int)mono_image_points.size();
	}

	bool CameraCalibrator::addImagePair(const std::string& left_image_path, const std::string& right_image_path)
	{
		std::vector<cv::Point2f> left_corners, right_corners;
		bool left_found = detectCorners(left_image_path, left_corners);
		bool right_found = detectCorners(right_image_path, right_corners);

		//DICe requires the target to be found in both cameras before an image set contributes
		//to the calibration (assemble_intersection_object_points() keeps only common points)
		if (!left_found || !right_found)
		{
			return false;
		}

		stereo_left_points.push_back(left_corners);
		stereo_right_points.push_back(right_corners);
		return true;
	}

	int CameraCalibrator::pairCount() const
	{
		return (int)stereo_left_points.size();
	}

	void CameraCalibrator::fillIntrinsics(const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs, CameraIntrinsics& intrinsics)
	{
		intrinsics.fx = (float)camera_matrix.at<double>(0, 0);
		intrinsics.fy = (float)camera_matrix.at<double>(1, 1);
		intrinsics.fs = 0.f; //cv::calibrateCamera does not solve for skew by default
		intrinsics.cx = (float)camera_matrix.at<double>(0, 2);
		intrinsics.cy = (float)camera_matrix.at<double>(1, 2);

		//OpenCorr's Calibration::distort() uses the same rational model and parameter
		//order as OpenCV: k1, k2, p1, p2, k3, k4, k5, k6 -- direct mapping, no reordering
		const double* dc = dist_coeffs.ptr<double>();
		size_t n = dist_coeffs.total();
		intrinsics.k1 = n > 0 ? (float)dc[0] : 0.f;
		intrinsics.k2 = n > 1 ? (float)dc[1] : 0.f;
		intrinsics.p1 = n > 2 ? (float)dc[2] : 0.f;
		intrinsics.p2 = n > 3 ? (float)dc[3] : 0.f;
		intrinsics.k3 = n > 4 ? (float)dc[4] : 0.f;
		intrinsics.k4 = n > 5 ? (float)dc[5] : 0.f;
		intrinsics.k5 = n > 6 ? (float)dc[6] : 0.f;
		intrinsics.k6 = n > 7 ? (float)dc[7] : 0.f;
	}

	float CameraCalibrator::calibrate(Calibration& camera)
	{
		if (mono_image_points.empty())
		{
			//addImage() silently excludes failed detections rather than treating them as an
			//error (matches DICe's own "exclude and continue" behavior) -- so it's reachable
			//in normal use for every image in a batch to fail (bad lighting, wrong
			//board_width/board_height). Without this check, cv::calibrateCamera's own
			//internal CV_Assert(nimages>0) throws an uncaught cv::Exception instead of a
			//clear, catchable error.
			throw std::string("CameraCalibrator::calibrate(): no successfully-detected images to calibrate from");
		}

		std::vector<std::vector<cv::Point3f>> object_points(mono_image_points.size(), object_points_template);

		cv::Mat camera_matrix, dist_coeffs;
		std::vector<cv::Mat> rvecs, tvecs;

		double rms = cv::calibrateCamera(object_points, mono_image_points, image_size,
			camera_matrix, dist_coeffs, rvecs, tvecs);

		camera.clear();
		fillIntrinsics(camera_matrix, dist_coeffs, camera.intrinsics);
		camera.updateMatrices();

		return (float)rms;
	}

	float CameraCalibrator::calibrateStereo(Calibration& left_camera, Calibration& right_camera, CameraExtrinsics& right_extrinsics)
	{
		if (stereo_left_points.empty())
		{
			throw std::string("CameraCalibrator::calibrateStereo(): no successfully-detected image pairs to calibrate from");
		}

		std::vector<std::vector<cv::Point3f>> object_points(stereo_left_points.size(), object_points_template);

		cv::Mat camera_matrix_l = cv::initCameraMatrix2D(object_points, stereo_left_points, image_size, 0);
		cv::Mat camera_matrix_r = cv::initCameraMatrix2D(object_points, stereo_right_points, image_size, 0);
		cv::Mat dist_l, dist_r, R, T, E, F;

		double rms = cv::stereoCalibrate(object_points, stereo_left_points, stereo_right_points,
			camera_matrix_l, dist_l, camera_matrix_r, dist_r, image_size, R, T, E, F,
			0, cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 1000, 1e-7));

		left_camera.clear();
		fillIntrinsics(camera_matrix_l, dist_l, left_camera.intrinsics);
		left_camera.updateMatrices(); //left camera stays at the identity/origin, right_extrinsics below is expressed relative to it

		right_camera.clear();
		fillIntrinsics(camera_matrix_r, dist_r, right_camera.intrinsics);

		cv::Mat rvec;
		cv::Rodrigues(R, rvec);
		right_extrinsics.rx = (float)rvec.at<double>(0);
		right_extrinsics.ry = (float)rvec.at<double>(1);
		right_extrinsics.rz = (float)rvec.at<double>(2);
		right_extrinsics.tx = (float)T.at<double>(0);
		right_extrinsics.ty = (float)T.at<double>(1);
		right_extrinsics.tz = (float)T.at<double>(2);
		right_camera.updateCalibration(right_camera.intrinsics, right_extrinsics);

		//calibration quality check, following DICe::Calibration::calibrate() (DICe_Calibration.cpp):
		//because the fundamental matrix implicitly encodes the full stereo geometry, the epipolar
		//constraint m2^T*F*m1=0 gives a per-point residual that is a genuinely different signal
		//than the RMS reprojection error above -- a low RMS can still hide a poorly-conditioned
		//region if the calibration images didn't cover it well
		epipolar_residuals.clear();
		epipolar_residuals.reserve(stereo_left_points.size());
		for (size_t i = 0; i < stereo_left_points.size(); i++)
		{
			cv::Mat pts_l(stereo_left_points[i]);
			cv::Mat pts_r(stereo_right_points[i]);
			cv::undistortPoints(pts_l, pts_l, camera_matrix_l, dist_l, cv::Mat(), camera_matrix_l);
			cv::undistortPoints(pts_r, pts_r, camera_matrix_r, dist_r, cv::Mat(), camera_matrix_r);

			std::vector<cv::Vec3f> lines_l, lines_r;
			cv::computeCorrespondEpilines(pts_l, 1, F, lines_l);
			cv::computeCorrespondEpilines(pts_r, 2, F, lines_r);

			double pair_error = 0.0;
			int npt = (int)stereo_left_points[i].size();
			for (int j = 0; j < npt; j++)
			{
				double err = std::fabs(stereo_left_points[i][j].x * lines_r[j][0] + stereo_left_points[i][j].y * lines_r[j][1] + lines_r[j][2])
					+ std::fabs(stereo_right_points[i][j].x * lines_l[j][0] + stereo_right_points[i][j].y * lines_l[j][1] + lines_l[j][2]);
				pair_error += err;
			}
			epipolar_residuals.push_back((float)(pair_error / npt));
		}

		return (float)rms;
	}

	const std::vector<float>& CameraCalibrator::epipolarResiduals() const
	{
		return epipolar_residuals;
	}

}//namespace opencorr
