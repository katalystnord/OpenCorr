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

#pragma once

#ifndef _CAMERA_CALIBRATOR_H_
#define _CAMERA_CALIBRATOR_H_

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>

#include "oc_calibration.h"

namespace opencorr
{
	//Solves for camera intrinsics/extrinsics from checkerboard target images
	//and populates OpenCorr's existing Calibration class -- OpenCorr's own
	//Calibration only consumes already-known parameters (matrix construction,
	//distort/undistort), it has no solve step at all. This wraps cv::
	//calibrateCamera/stereoCalibrate (both already reachable through
	//OpenCorr's existing OpenCV dependency) plus the stereo calibration-
	//quality metric from DICe::Calibration::calibrate() (dicengine/dice,
	//DICe_Calibration.cpp, BSD-3-Clause): the epipolar-constraint residual
	//m2^T*F*m1 averaged per image pair, which is a genuinely different
	//failure-mode signal than the overall RMS reprojection error -- a pair
	//can have a low RMS while still being poorly conditioned for triangulation
	//in one region of the image, and the per-pair epipolar residual surfaces
	//that where a single aggregate RMS number does not.
	//
	//Scope note: this covers the checkerboard target type only. DICe also
	//supports dot targets with "donut" marker dots for axis/origin
	//determination (SimpleBlobDetector + a nanoflann grid-assembly step) --
	//a substantially larger, separate piece of work, intentionally not
	//included here.

	class CameraCalibrator
	{
	public:
		CameraCalibrator(int board_width, int board_height, float square_size);
		~CameraCalibrator();

		//detect checkerboard corners in a single image; returns true on success.
		//failed images are simply not added (not counted as an error) -- mirrors
		//DICe's own "exclude the image and continue" behavior in extract_checkerboard_intersections()
		bool addImage(const std::string& image_path);
		int imageCount() const;

		//run single-camera calibration on all successfully-added images;
		//fills camera.intrinsics and updates its matrices; returns the RMS reprojection error
		float calibrate(Calibration& camera);

		//detect checkerboard corners in a left/right image pair; returns true only if found in both
		bool addImagePair(const std::string& left_image_path, const std::string& right_image_path);
		int pairCount() const;

		//run stereo calibration on all successfully-added pairs; fills left/right camera
		//intrinsics (left_camera's extrinsics are left at the identity/origin -- the world
		//origin is defined at the left camera, matching the convention right_extrinsics
		//is expressed relative to) and right_extrinsics (right camera's pose relative to
		//the left camera); also computes the per-pair epipolar residuals (see below);
		//returns the overall RMS reprojection error
		float calibrateStereo(Calibration& left_camera, Calibration& right_camera, CameraExtrinsics& right_extrinsics);

		//one entry per successfully-added stereo pair (same order as addImagePair calls),
		//average |m2^T F m1| over that pair's corner correspondences -- only populated
		//after calibrateStereo() has been called
		const std::vector<float>& epipolarResiduals() const;

	private:
		int board_width, board_height; //number of internal corners along each dimension
		float square_size; //physical size of one checkerboard square, in whatever unit the caller wants results expressed in
		std::vector<cv::Point3f> object_points_template; //3D checkerboard corners in target coordinates (z=0), built once from board_width/height/square_size

		cv::Size image_size;
		bool image_size_set;

		std::vector<std::vector<cv::Point2f>> mono_image_points;
		std::vector<std::vector<cv::Point2f>> stereo_left_points;
		std::vector<std::vector<cv::Point2f>> stereo_right_points;

		std::vector<float> epipolar_residuals;

		bool detectCorners(const std::string& image_path, std::vector<cv::Point2f>& corners);
		static void fillIntrinsics(const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs, CameraIntrinsics& intrinsics);
	};

}//namespace opencorr

#endif //_CAMERA_CALIBRATOR_H_
