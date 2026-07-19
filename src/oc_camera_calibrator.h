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

/*
 * Portions of this file are ported from DICe (dicengine/dice), BSD-3-Clause.
 * See THIRD-PARTY-LICENSES.md for the required copyright notice, list of
 * conditions, and disclaimer.
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
	//Scope note: covers checkerboard targets (cv::findChessboardCorners) and
	//dot targets with 3 "donut" marker dots for axis/origin determination
	//(DICe::Calibration::extract_dot_target_points(), DICe_OpenCVServerUtils.cpp
	//opencv_dot_targets()/get_dot_markers()/reorder_keypoints() -- dicengine/dice,
	//BSD-3-Clause). Unlike checkerboard detection (all-or-nothing: every internal
	//corner or none), dot-target detection is inherently partial -- individual
	//dots can be missed to lighting/occlusion -- so each detected image/pair
	//carries its own subset of object_points_template rather than assuming the
	//full grid, mirroring DICe's own per-image "common point" handling.

	enum class CalibrationTargetType
	{
		CHECKER_BOARD,
		DOT_TARGET //black dots on white background, 3 donut markers mark origin/x-axis/y-axis
	};

	class CameraCalibrator
	{
	public:
		CameraCalibrator(int board_width, int board_height, float square_size);

		//dot-target constructor. origin_x/origin_y is the grid index (0-based) of the
		//origin marker dot; num_fiducials_origin_to_x/y_marker is the number of grid
		//points (inclusive) from the origin marker to the x-axis/y-axis marker dot along
		//each axis -- mirrors DICe's origin_loc_x_/origin_loc_y_/
		//num_fiducials_origin_to_x_marker_/num_fiducials_origin_to_y_marker_
		CameraCalibrator(int num_fiducials_x, int num_fiducials_y, float dot_spacing,
			int origin_x, int origin_y, int num_fiducials_origin_to_x_marker, int num_fiducials_origin_to_y_marker);

		~CameraCalibrator();

		//detect the calibration target in a single image; returns true on success.
		//failed images are simply not added (not counted as an error) -- mirrors
		//DICe's own "exclude the image and continue" behavior in extract_checkerboard_intersections()
		//and extract_dot_target_points()
		bool addImage(const std::string& image_path);
		int imageCount() const;

		//run single-camera calibration on all successfully-added images;
		//fills camera.intrinsics and updates its matrices; returns the RMS reprojection error
		float calibrate(Calibration& camera);

		//detect the calibration target in a left/right image pair; returns true only if
		//found in both, and (for dot targets) at least one point is common to both images
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
		CalibrationTargetType target_type = CalibrationTargetType::CHECKER_BOARD;

		int board_width, board_height; //checkerboard: internal corners along each dimension; dot target: fiducials along each dimension
		float square_size; //checkerboard: square size; dot target: dot spacing -- physical units, whatever unit the caller wants results expressed in
		std::vector<cv::Point3f> object_points_template; //3D target points in target coordinates (z=0), row-major (index = row*board_width + col)

		//dot-target-only parameters, mirroring DICe's origin_loc_x_/origin_loc_y_/
		//num_fiducials_origin_to_x_marker_/num_fiducials_origin_to_y_marker_
		int dot_origin_x = 0, dot_origin_y = 0;
		int dot_origin_to_x_marker = 0, dot_origin_to_y_marker = 0;

		cv::Size image_size;
		bool image_size_set;

		//each image/pair keeps its own subset of object_points_template (grid_indices
		//identify which template entries a detection's points correspond to) -- for
		//checkerboard this subset is always the full grid in order, for dot targets
		//it is whatever subset was actually detected
		std::vector<std::vector<cv::Point2f>> mono_image_points;
		std::vector<std::vector<int>> mono_grid_indices;
		std::vector<std::vector<cv::Point2f>> stereo_left_points;
		std::vector<std::vector<cv::Point2f>> stereo_right_points;
		std::vector<std::vector<int>> stereo_grid_indices;

		std::vector<float> epipolar_residuals;

		bool detectCorners(const std::string& image_path, std::vector<cv::Point2f>& corners);

		//detects a dot target; dots[k] is the image position of the grid point at
		//flattened index grid_indices[k] (= row*board_width + col) in object_points_template,
		//including the 3 donut marker dots at their known grid indices. Returns false
		//(without throwing) if the target, or the 3 donut markers specifically, aren't found.
		bool detectDots(const std::string& image_path, std::vector<cv::Point2f>& dots, std::vector<int>& grid_indices);

		//loads image_path (grayscale) and checks/sets image_size; shared by detectCorners/detectDots
		bool loadAndCheckSize(const std::string& image_path, cv::Mat& img);

		static void fillIntrinsics(const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs, CameraIntrinsics& intrinsics);

		//builds the object-point subset for a detection given its grid_indices
		std::vector<cv::Point3f> objectPointsFor(const std::vector<int>& grid_indices) const;
	};

}//namespace opencorr

#endif //_CAMERA_CALIBRATOR_H_
