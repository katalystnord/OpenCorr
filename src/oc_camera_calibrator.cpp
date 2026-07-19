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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

#include "oc_camera_calibrator.h"
#include "oc_nearest_neighbor.h"

namespace opencorr
{
	namespace
	{
		//builds a blob detector for the given minimum blob area. Detector construction only
		//depends on min_blob_size, which stays constant across an entire threshold-sweep
		//attempt in detectDots() -- callers build one of these once per attempt and reuse it
		//across all threshold iterations instead of rebuilding per call.
		//
		//filterByColor stays at its default (true, blobColor=0 -- select dark blobs specifically):
		//this is load-bearing, not incidental -- detectDots() relies on the SAME threshold
		//producing a dark-blob-only search on not_src (finds donut holes, ~3 results) and on
		//bi_src (finds plain dots, ~num_grid results); disabling color filtering would make both
		//passes match bright and dark blobs alike, breaking that separation.
		//
		//filterByInertia/filterByConvexity ARE disabled: their strict defaults (minInertiaRatio
		//~0.1, minConvexity ~0.95) are shape-strictness knobs tuned for near-perfect circles and
		//can reject legitimate dots on real (noisy/compressed/imperfectly printed) camera images
		//-- DICe's own equivalent detector setup avoids exactly this by leaving these off. This
		//doesn't affect which blobs get selected by brightness, only how strict their shape has
		//to be, so it's safe to relax independently of filterByColor.
		cv::Ptr<cv::SimpleBlobDetector> makeBlobDetector(int min_blob_size)
		{
			cv::SimpleBlobDetector::Params params;
			params.filterByArea = true;
			params.minArea = (float)min_blob_size;
			params.maxArea = 10e4f;
			params.filterByInertia = false;
			params.filterByConvexity = false;
			return cv::SimpleBlobDetector::create(params);
		}

		//detects blob keypoints in img at the given threshold. When donut_test is true,
		//keeps only blobs whose connected-component "shadow" splits into exactly 3 zones
		//along a horizontal scan through the blob center (background/ring/hole) --
		//identifying the "donut" marker dots vs. plain filled dots. invert selects which of
		//the thresholded image / its complement the blob detector runs on (true: dots are
		//darker than the background). Mirrors DICe's get_dot_markers() (DICe_OpenCVServerUtils.cpp).
		std::vector<cv::KeyPoint> findDotMarkers(const cv::Mat& img, int thresh, bool invert, bool donut_test,
			const cv::Ptr<cv::SimpleBlobDetector>& detector)
		{
			cv::Mat bi_src, not_src;
			cv::threshold(img, bi_src, thresh, 255, cv::THRESH_BINARY);
			cv::bitwise_not(bi_src, not_src);

			std::vector<cv::KeyPoint> keypoints;
			if (invert)
			{
				detector->detect(not_src, keypoints);
			}
			else
			{
				detector->detect(bi_src, keypoints);
			}

			if (keypoints.empty())
			{
				return keypoints;
			}

			//connected-component labeling is only needed for the donut zone-count test below --
			//computed here (not unconditionally at the top of the function) so the general
			//(non-donut) dot-detection pass, and every threshold-sweep iteration that finds no
			//keypoints at all, don't pay for a full-image labeling pass whose result would be
			//immediately discarded
			cv::Mat label_image;
			cv::Mat stats, centroids;
			int num_labels = 0;
			if (donut_test)
			{
				if (invert)
				{
					num_labels = cv::connectedComponentsWithStats(bi_src, label_image, stats, centroids, 8, CV_32S);
				}
				else
				{
					num_labels = cv::connectedComponentsWithStats(not_src, label_image, stats, centroids, 8, CV_32S);
				}
			}

			float avg_size = 0.f;
			for (auto& kp : keypoints)
			{
				avg_size += kp.size;
			}
			avg_size /= (float)keypoints.size();

			if (!donut_test)
			{
				for (int i = (int)keypoints.size() - 1; i >= 0; i--)
				{
					if (std::fabs(keypoints[i].size - avg_size) / avg_size > 0.5f)
					{
						keypoints.erase(keypoints.begin() + i);
					}
				}
				return keypoints;
			}

			const int donut_span = std::max(1, (int)(2.f * avg_size));
			const int width = img.cols;
			for (int i = (int)keypoints.size() - 1; i >= 0; i--)
			{
				int cx = (int)keypoints[i].pt.x;
				int cy = (int)keypoints[i].pt.y;
				int x_lo = cx - donut_span;
				int x_hi = cx + donut_span;
				if (x_lo < 0 || x_hi >= width || cy < 0 || cy >= img.rows ||
					label_image.at<int>(cy, x_lo) != label_image.at<int>(cy, x_hi))
				{
					keypoints.erase(keypoints.begin() + i);
					continue;
				}
				std::set<int> zones;
				for (int x = x_lo; x < x_hi; x++)
				{
					zones.insert(label_image.at<int>(cy, x));
				}
				if (zones.size() != 3)
				{
					keypoints.erase(keypoints.begin() + i);
				}
			}

			//re-derive keypoint size from its connected-component area and drop outliers a
			//second time. keypoints[i].pt sits at the detected hole's centroid (that's what
			//makes it a donut detection), and label_image's positive labels here come from the
			//hole (an isolated component, disconnected from the outer background by the ring)
			//-- so this measures the HOLE's precise pixel-count area, not the whole ring+hole
			//marker's area. That's still a real improvement over SimpleBlobDetector's own
			//kp.size (a coarse diameter estimate of the same hole): a true area comparison
			//across the 3 markers' holes is a more precise per-marker consistency check,
			//assuming the 3 markers are printed with matching hole sizes
			if (keypoints.size() == 3 && num_labels >= 3)
			{
				std::vector<int> areas(keypoints.size());
				float avg_area = 0.f;
				for (size_t i = 0; i < keypoints.size(); i++)
				{
					int label = label_image.at<int>((int)keypoints[i].pt.y, (int)keypoints[i].pt.x);
					areas[i] = stats.at<int>(label, cv::CC_STAT_AREA);
					avg_area += (float)areas[i];
				}
				avg_area /= (float)keypoints.size();
				for (int i = (int)keypoints.size() - 1; i >= 0; i--)
				{
					if (areas[i] <= 0 || std::fabs((float)areas[i] - avg_area) / avg_area > 2.f)
					{
						keypoints.erase(keypoints.begin() + i);
					}
				}
			}

			return keypoints;
		}

		//orders 3 raw donut keypoints into [origin, x-axis marker, y-axis marker]. For each
		//of the 3 lines through a pair of markers, counts how many detected dots lie close to
		//it -- the two longest fiducial rows/columns produce the highest counts. The origin is
		//the marker common to the two highest-count lines; the remaining two markers are the
		//x-axis end (on the higher-count of the two remaining lines) and y-axis end.
		//Mirrors DICe's reorder_keypoints() (DICe_OpenCVServerUtils.cpp).
		bool reorderMarkers(std::vector<cv::KeyPoint>& markers, const std::vector<cv::KeyPoint>& dots)
		{
			//line[i] passes through the two markers NOT at index i. Point-to-line distance via
			//the 2D cross product |cross(p2-p1, q-p1)| / |p2-p1| avoids the slope-based
			//(y2-y1)/(x2-x1) formula's singularity: two markers sharing (or nearly sharing) an
			//image x-coordinate -- a real, non-exotic case for a near-fronto-parallel
			//calibration pose with a near-vertical fiducial column -- would otherwise produce an
			//Inf/NaN slope, silently zeroing that line's vote instead of the intended failure path
			cv::Point2f line_p1[3] = { markers[1].pt, markers[0].pt, markers[0].pt };
			cv::Point2f line_p2[3] = { markers[2].pt, markers[2].pt, markers[1].pt };
			float count[3] = { 0.f, 0.f, 0.f };

			const float dist_tol = 5.f;
			for (int i = 0; i < 3; i++)
			{
				cv::Point2f dir = line_p2[i] - line_p1[i];
				float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
				if (len < 1e-3f)
				{
					continue; //degenerate (coincident) marker pair -- leave this line's count at 0
				}
				for (auto& d : dots)
				{
					cv::Point2f rel = d.pt - line_p1[i];
					float cross = dir.x * rel.y - dir.y * rel.x;
					float dist = std::fabs(cross) / len;
					if (dist < dist_tol)
					{
						count[i] += 1.f;
					}
				}
			}

			int order[3] = { 0, 1, 2 };
			std::sort(order, order + 3, [&](int a, int b2) { return count[a] > count[b2]; });

			std::set<int> side[3] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };

			int origin_idx = -1;
			for (int a : side[order[0]])
			{
				for (int b2 : side[order[1]])
				{
					if (a == b2)
					{
						origin_idx = a;
					}
				}
			}
			if (origin_idx < 0)
			{
				return false; //degenerate marker geometry (e.g. collinear markers)
			}

			side[order[0]].erase(origin_idx);
			side[order[1]].erase(origin_idx);
			if (side[order[0]].empty() || side[order[1]].empty())
			{
				return false;
			}
			int x_idx = *side[order[0]].begin();
			int y_idx = *side[order[1]].begin();
			if (x_idx == y_idx)
			{
				return false;
			}

			std::vector<cv::KeyPoint> ordered = { markers[origin_idx], markers[x_idx], markers[y_idx] };
			markers = ordered;
			return true;
		}
	}

	CameraCalibrator::CameraCalibrator(int board_width, int board_height, float square_size)
		: target_type(CalibrationTargetType::CHECKER_BOARD),
		board_width(board_width), board_height(board_height), square_size(square_size), image_size_set(false)
	{
		for (int r = 0; r < board_height; r++)
		{
			for (int c = 0; c < board_width; c++)
			{
				object_points_template.push_back(cv::Point3f(c * square_size, r * square_size, 0.f));
			}
		}
	}

	CameraCalibrator::CameraCalibrator(int num_fiducials_x, int num_fiducials_y, float dot_spacing,
		int origin_x, int origin_y, int num_fiducials_origin_to_x_marker, int num_fiducials_origin_to_y_marker)
		: target_type(CalibrationTargetType::DOT_TARGET),
		board_width(num_fiducials_x), board_height(num_fiducials_y), square_size(dot_spacing),
		dot_origin_x(origin_x), dot_origin_y(origin_y),
		dot_origin_to_x_marker(num_fiducials_origin_to_x_marker), dot_origin_to_y_marker(num_fiducials_origin_to_y_marker),
		image_size_set(false)
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

	bool CameraCalibrator::loadAndCheckSize(const std::string& image_path, cv::Mat& img)
	{
		img = cv::imread(image_path, cv::IMREAD_GRAYSCALE);
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

		return true;
	}

	bool CameraCalibrator::detectCorners(const std::string& image_path, std::vector<cv::Point2f>& corners)
	{
		cv::Mat img;
		if (!loadAndCheckSize(image_path, img))
		{
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

	bool CameraCalibrator::detectDots(const std::string& image_path, std::vector<cv::Point2f>& dots_out, std::vector<int>& grid_indices_out)
	{
		dots_out.clear();
		grid_indices_out.clear();

		cv::Mat img;
		if (!loadAndCheckSize(image_path, img))
		{
			return false;
		}

		const int num_fiducials_x = board_width;
		const int num_fiducials_y = board_height;
		const int num_grid = num_fiducials_x * num_fiducials_y;

		//--- phase 1: sweep thresholds to find exactly 3 donut markers, retrying once at a
		//    smaller min_blob_size if the first pass fails entirely ---
		const int threshold_start = 20, threshold_end = 250, threshold_step = 5;
		int min_blob_size = 100;
		std::vector<cv::KeyPoint> markers;
		bool found_markers = false;
		int good_thresh = threshold_start;

		for (int attempt = 0; attempt < 2 && !found_markers; attempt++)
		{
			if (attempt == 1)
			{
				min_blob_size = 10;
			}

			//min_blob_size (the only detector parameter that varies) is constant across this
			//whole attempt's threshold sweep, so the detector is built once per attempt
			//instead of once per threshold iteration
			cv::Ptr<cv::SimpleBlobDetector> detector = makeBlobDetector(min_blob_size);

			int thresh_first = 0, thresh_last = 0;
			for (int t = threshold_start; t <= threshold_end; t += threshold_step)
			{
				markers = findDotMarkers(img, t, true, true, detector);
				if (markers.size() == 3)
				{
					if (thresh_first == 0)
					{
						thresh_first = t;
					}
					thresh_last = t;
				}
				else if (thresh_first != 0)
				{
					break;
				}
			}

			if (thresh_first != 0)
			{
				good_thresh = (thresh_first + thresh_last) / 2;
				good_thresh -= good_thresh % threshold_step;
				markers = findDotMarkers(img, good_thresh, true, true, detector);
				found_markers = (markers.size() == 3);
			}
		}

		if (!found_markers)
		{
			std::cerr << "CameraCalibrator: could not identify 3 donut marker dots, excluding: " << image_path << std::endl;
			return false;
		}

		//--- phase 2: determine a working threshold from the gray level between two markers,
		//    then detect all dots at that threshold ---
		float du = markers[1].pt.x - markers[0].pt.x;
		float dv = markers[1].pt.y - markers[0].pt.y;
		float mag = std::sqrt(du * du + dv * dv);
		if (mag < 1.f)
		{
			std::cerr << "CameraCalibrator: degenerate marker geometry, excluding: " << image_path << std::endl;
			return false;
		}
		float nu = du / mag, nv = dv / mag;

		int min_gray = 255, max_gray = 0;
		for (int r = 1; r < (int)mag; r++)
		{
			int ix = (int)(markers[0].pt.x + r * nu);
			int iy = (int)(markers[0].pt.y + r * nv);
			if (ix < 0 || ix >= img.cols || iy < 0 || iy >= img.rows)
			{
				continue;
			}
			int g = img.at<uchar>(iy, ix);
			min_gray = std::min(min_gray, g);
			max_gray = std::max(max_gray, g);
		}
		int dots_thresh = (min_gray + max_gray) / 2;

		//the general dot pass uses the OPPOSITE invert setting from the marker pass: markers
		//are identified by their hole showing up as an isolated blob on the inverted image
		//(not_src), while plain dots (solid, no hole) are found directly on the un-inverted
		//binary image (bi_src) -- confirmed empirically (see dot_target_smoke_test.cpp), also
		//matches DICe's own get_dot_markers(img_cpy, dots, i_thresh, !invert, ...) call
		cv::Ptr<cv::SimpleBlobDetector> general_detector = makeBlobDetector(min_blob_size);
		std::vector<cv::KeyPoint> dots = findDotMarkers(img, dots_thresh, false, false, general_detector);
		if ((int)dots.size() < (int)(num_grid * 0.7))
		{
			//the between-markers gray level didn't produce enough dots -- fall back to the
			//threshold that worked for marker detection
			dots = findDotMarkers(img, good_thresh, false, false, general_detector);
		}
		if (dots.empty())
		{
			std::cerr << "CameraCalibrator: zero dots found, excluding: " << image_path << std::endl;
			return false;
		}

		//--- phase 3: identify origin/x-marker/y-marker among the 3 raw markers ---
		if (!reorderMarkers(markers, dots))
		{
			std::cerr << "CameraCalibrator: could not order marker dots into origin/x-axis/y-axis, excluding: " << image_path << std::endl;
			return false;
		}

		//--- phase 4: nearest-neighbor grid walk outward from the markers, using OpenCorr's
		//    own NearestNeighbor (nanoflann-backed) instead of a hand-rolled kd-tree ---
		std::vector<Point2D> cloud_points;
		cloud_points.reserve(dots.size() + 3);
		for (auto& d : dots)
		{
			cloud_points.push_back(Point2D(d.pt.x, d.pt.y));
		}
		for (auto& mk : markers)
		{
			cloud_points.push_back(Point2D(mk.pt.x, mk.pt.y));
		}

		NearestNeighbor nn;
		nn.assignPoints(cloud_points);
		nn.setSearchK(1);
		nn.constructKdTree();

		std::vector<int> use_count(cloud_points.size(), 0);
		std::vector<cv::Point2f> dot_grid(num_grid, cv::Point2f(-1.f, -1.f));

		if (dot_origin_to_x_marker <= 1 || dot_origin_to_y_marker <= 1)
		{
			std::cerr << "CameraCalibrator: num_fiducials_origin_to_x/y_marker must be > 1" << std::endl;
			return false;
		}
		if (dot_origin_x < 0 || dot_origin_x >= num_fiducials_x || dot_origin_y < 0 || dot_origin_y >= num_fiducials_y)
		{
			std::cerr << "CameraCalibrator: origin_x/origin_y must be within [0, num_fiducials_x/y)" << std::endl;
			return false;
		}
		float udx = (markers[1].pt.x - markers[0].pt.x) / (float)(dot_origin_to_x_marker - 1);
		float udy = (markers[1].pt.y - markers[0].pt.y) / (float)(dot_origin_to_x_marker - 1);
		float vdx = (markers[2].pt.x - markers[0].pt.x) / (float)(dot_origin_to_y_marker - 1);
		float vdy = (markers[2].pt.y - markers[0].pt.y) / (float)(dot_origin_to_y_marker - 1);
		//separate acceptance radii for the two axes -- a single shared threshold derived from
		//only one axis's spacing is wrong whenever pitch/foreshortening differs between the
		//origin-to-x-marker and origin-to-y-marker directions (non-square dot spacing,
		//asymmetric marker placement, or ordinary perspective tilt)
		float dist_threshold_u = 0.25f * udx * udx + 0.25f * udy * udy;
		float dist_threshold_v = 0.25f * vdx * vdx + 0.25f * vdy * vdy;

		auto tryPlace = [&](int gi, float qx, float qy, bool boundary_check, float dist_threshold) -> bool
		{
			if (gi < 0 || gi >= num_grid || dot_grid[gi].x >= 0.f)
			{
				return false;
			}
			if (boundary_check && (qx < 20.f || qx > image_size.width - 20.f))
			{
				return false;
			}
			std::vector<uint32_t> idx;
			std::vector<float> dist2;
			if (nn.knnSearch(Point3D(qx, qy, 0.f), 1, idx, dist2) < 1)
			{
				return false;
			}
			if (dist2[0] >= dist_threshold || use_count[idx[0]] > 0)
			{
				return false;
			}
			dot_grid[gi] = cv::Point2f(cloud_points[idx[0]].x, cloud_points[idx[0]].y);
			use_count[idx[0]]++;
			return true;
		};

		//x-axis, at/before the origin
		for (int i = 0; i <= dot_origin_x; i++)
		{
			float qx = (i - dot_origin_x) * udx + markers[0].pt.x;
			float qy = (i - dot_origin_x) * udy + markers[0].pt.y;
			tryPlace(dot_origin_y * num_fiducials_x + i, qx, qy, true, dist_threshold_u);
		}
		if (dot_grid[dot_origin_y * num_fiducials_x + dot_origin_x].x < 0.f)
		{
			std::cerr << "CameraCalibrator: could not place origin dot, excluding: " << image_path << std::endl;
			return false;
		}
		//x-axis, past the origin, extrapolating from the last two placed points
		for (int i = dot_origin_x + 1; i < num_fiducials_x; i++)
		{
			int prev = dot_origin_y * num_fiducials_x + i - 1;
			if (dot_grid[prev].x < 0.f)
			{
				break;
			}
			float qx, qy;
			if (i == dot_origin_x + 1)
			{
				qx = dot_grid[prev].x + udx;
				qy = dot_grid[prev].y + udy;
			}
			else
			{
				int prev2 = dot_origin_y * num_fiducials_x + i - 2;
				qx = 2.f * dot_grid[prev].x - dot_grid[prev2].x;
				qy = 2.f * dot_grid[prev].y - dot_grid[prev2].y;
			}
			tryPlace(dot_origin_y * num_fiducials_x + i, qx, qy, false, dist_threshold_u);
		}
		//one row up from the x-axis
		for (int i = 0; i < num_fiducials_x; i++)
		{
			int base = dot_origin_y * num_fiducials_x + i;
			if (dot_grid[base].x < 0.f)
			{
				continue;
			}
			tryPlace((dot_origin_y + 1) * num_fiducials_x + i, dot_grid[base].x + vdx, dot_grid[base].y + vdy, false, dist_threshold_v);
		}
		//remaining rows above, extrapolating from the two rows below
		for (int j = dot_origin_y + 2; j < num_fiducials_y; j++)
		{
			for (int i = 0; i < num_fiducials_x; i++)
			{
				int r1 = (j - 1) * num_fiducials_x + i, r2 = (j - 2) * num_fiducials_x + i;
				if (dot_grid[r1].x < 0.f || dot_grid[r2].x < 0.f)
				{
					continue;
				}
				tryPlace(j * num_fiducials_x + i, 2.f * dot_grid[r1].x - dot_grid[r2].x, 2.f * dot_grid[r1].y - dot_grid[r2].y, false, dist_threshold_v);
			}
		}
		//rows below the x-axis, extrapolating downward
		for (int j = dot_origin_y - 1; j >= 0; j--)
		{
			for (int i = 0; i < num_fiducials_x; i++)
			{
				int r1 = (j + 1) * num_fiducials_x + i, r2 = (j + 2) * num_fiducials_x + i;
				if (dot_grid[r1].x < 0.f || dot_grid[r2].x < 0.f)
				{
					continue;
				}
				tryPlace(j * num_fiducials_x + i, 2.f * dot_grid[r1].x - dot_grid[r2].x, 2.f * dot_grid[r1].y - dot_grid[r2].y, false, dist_threshold_v);
			}
		}

		//--- phase 5: stamp the 3 markers at their known grid indices (exact positions,
		//    overriding whatever the grid walk may have placed there), assemble output ---
		int marker_x_col = dot_origin_x + dot_origin_to_x_marker - 1;
		int marker_y_row = dot_origin_y + dot_origin_to_y_marker - 1;
		//validate the ROW/COLUMN components individually, not just the flattened total -- a
		//column overflow (marker_x_col >= num_fiducials_x) can still flatten to a value inside
		//[0, num_grid) by wrapping into the next row, silently landing on and overwriting an
		//unrelated, already-placed grid cell instead of failing
		if (marker_x_col < 0 || marker_x_col >= num_fiducials_x || marker_y_row < 0 || marker_y_row >= num_fiducials_y)
		{
			std::cerr << "CameraCalibrator: marker grid position out of range, check num_fiducials_origin_to_x/y_marker" << std::endl;
			return false;
		}
		int marker_origin_idx = dot_origin_y * num_fiducials_x + dot_origin_x;
		int marker_x_idx = dot_origin_y * num_fiducials_x + marker_x_col;
		int marker_y_idx = marker_y_row * num_fiducials_x + dot_origin_x;
		dot_grid[marker_origin_idx] = markers[0].pt;
		dot_grid[marker_x_idx] = markers[1].pt;
		dot_grid[marker_y_idx] = markers[2].pt;

		for (int gi = 0; gi < num_grid; gi++)
		{
			if (dot_grid[gi].x < 0.f)
			{
				continue;
			}
			dots_out.push_back(dot_grid[gi]);
			grid_indices_out.push_back(gi);
		}

		if ((int)dots_out.size() < (int)(num_grid * 0.7))
		{
			std::cerr << "CameraCalibrator: not enough dots found (" << dots_out.size() << "/" << num_grid
				<< "), excluding: " << image_path << std::endl;
			dots_out.clear();
			grid_indices_out.clear();
			return false;
		}

		return true;
	}

	std::vector<cv::Point3f> CameraCalibrator::objectPointsFor(const std::vector<int>& grid_indices) const
	{
		std::vector<cv::Point3f> pts;
		pts.reserve(grid_indices.size());
		for (int gi : grid_indices)
		{
			pts.push_back(object_points_template[gi]);
		}
		return pts;
	}

	bool CameraCalibrator::addImage(const std::string& image_path)
	{
		std::vector<cv::Point2f> points;
		std::vector<int> grid_indices;
		bool found;

		if (target_type == CalibrationTargetType::DOT_TARGET)
		{
			found = detectDots(image_path, points, grid_indices);
		}
		else
		{
			found = detectCorners(image_path, points);
			if (found)
			{
				grid_indices.resize(object_points_template.size());
				std::iota(grid_indices.begin(), grid_indices.end(), 0);
			}
		}

		if (!found)
		{
			return false;
		}

		mono_image_points.push_back(points);
		mono_grid_indices.push_back(grid_indices);
		return true;
	}

	int CameraCalibrator::imageCount() const
	{
		return (int)mono_image_points.size();
	}

	bool CameraCalibrator::addImagePair(const std::string& left_image_path, const std::string& right_image_path)
	{
		std::vector<cv::Point2f> left_points, right_points;
		std::vector<int> left_indices, right_indices;
		bool left_found, right_found;

		if (target_type == CalibrationTargetType::DOT_TARGET)
		{
			left_found = detectDots(left_image_path, left_points, left_indices);
			right_found = detectDots(right_image_path, right_points, right_indices);
		}
		else
		{
			left_found = detectCorners(left_image_path, left_points);
			right_found = detectCorners(right_image_path, right_points);
			if (left_found)
			{
				left_indices.resize(object_points_template.size());
				std::iota(left_indices.begin(), left_indices.end(), 0);
			}
			right_indices = left_indices;
		}

		//DICe requires the target to be found in both cameras before an image set contributes
		//to the calibration (assemble_intersection_object_points() keeps only common points)
		if (!left_found || !right_found)
		{
			return false;
		}

		std::vector<cv::Point2f> common_left, common_right;
		std::vector<int> common_indices;
		for (size_t i = 0; i < left_indices.size(); i++)
		{
			for (size_t j = 0; j < right_indices.size(); j++)
			{
				if (left_indices[i] == right_indices[j])
				{
					common_left.push_back(left_points[i]);
					common_right.push_back(right_points[j]);
					common_indices.push_back(left_indices[i]);
					break;
				}
			}
		}

		//require the same 70%-of-grid floor detectDots() applies per image (mirrors DICe's own
		//per-image-set common-point tolerance in extract_dot_target_points()) -- a nonzero-only
		//floor isn't enough: two independently-thresholded detections can each individually
		//clear 70% while sharing very few common points (differing lighting/occlusion per
		//camera), and cv::stereoCalibrate/cv::initCameraMatrix2D need enough points per view for
		//a stable pose solve, not just a nonempty one
		if ((int)common_indices.size() < (int)(object_points_template.size() * 0.7))
		{
			std::cerr << "CameraCalibrator: too few common points between left/right images (" << common_indices.size()
				<< "/" << object_points_template.size() << "), excluding pair: "
				<< left_image_path << " / " << right_image_path << std::endl;
			return false;
		}

		stereo_left_points.push_back(common_left);
		stereo_right_points.push_back(common_right);
		stereo_grid_indices.push_back(common_indices);
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

		std::vector<std::vector<cv::Point3f>> object_points;
		object_points.reserve(mono_image_points.size());
		for (auto& indices : mono_grid_indices)
		{
			object_points.push_back(objectPointsFor(indices));
		}

		cv::Mat camera_matrix, dist_coeffs;
		std::vector<cv::Mat> rvecs, tvecs;

		double rms = cv::calibrateCamera(object_points, mono_image_points, image_size,
			camera_matrix, dist_coeffs, rvecs, tvecs);

		camera.clear();
		fillIntrinsics(camera_matrix, dist_coeffs, camera.intrinsics);
		camera.updateMatrices();
		//every existing caller of Calibration in this codebase (Stereovision, EpipolarSearch,
		//and every example that builds one) calls prepare() immediately after loading
		//intrinsics/extrinsics, before the Calibration is used -- it populates the
		//undistort() maps that Stereovision::reconstruct() indexes. Without this, a
		//Calibration produced here and handed directly to Stereovision/EpipolarSearch would
		//index an empty (0x0) map.
		camera.prepare(image_size.height, image_size.width);

		return (float)rms;
	}

	float CameraCalibrator::calibrateStereo(Calibration& left_camera, Calibration& right_camera, CameraExtrinsics& right_extrinsics)
	{
		if (stereo_left_points.empty())
		{
			throw std::string("CameraCalibrator::calibrateStereo(): no successfully-detected image pairs to calibrate from");
		}

		std::vector<std::vector<cv::Point3f>> object_points;
		object_points.reserve(stereo_left_points.size());
		for (auto& indices : stereo_grid_indices)
		{
			object_points.push_back(objectPointsFor(indices));
		}

		cv::Mat camera_matrix_l = cv::initCameraMatrix2D(object_points, stereo_left_points, image_size, 0);
		cv::Mat camera_matrix_r = cv::initCameraMatrix2D(object_points, stereo_right_points, image_size, 0);
		cv::Mat dist_l, dist_r, R, T, E, F;

		double rms = cv::stereoCalibrate(object_points, stereo_left_points, stereo_right_points,
			camera_matrix_l, dist_l, camera_matrix_r, dist_r, image_size, R, T, E, F,
			0, cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 1000, 1e-7));

		left_camera.clear();
		fillIntrinsics(camera_matrix_l, dist_l, left_camera.intrinsics);
		left_camera.updateMatrices(); //left camera stays at the identity/origin, right_extrinsics below is expressed relative to it
		left_camera.prepare(image_size.height, image_size.width); //see calibrate()'s own comment on why this is required

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
		right_camera.prepare(image_size.height, image_size.width);

		//calibration quality check, following DICe::Calibration::calibrate() (DICe_Calibration.cpp):
		//because the fundamental matrix implicitly encodes the full stereo geometry, the epipolar
		//constraint m2^T*F*m1=0 gives a per-point residual that is a genuinely different signal
		//than the RMS reprojection error above -- a low RMS can still hide a poorly-conditioned
		//region if the calibration images didn't cover it well
		epipolar_residuals.clear();
		epipolar_residuals.reserve(stereo_left_points.size());
		for (size_t i = 0; i < stereo_left_points.size(); i++)
		{
			//undistort into NEW vectors rather than in place: cv::Mat's constructor from a
			//std::vector<Point2f> wraps that vector's own buffer without copying (verified),
			//so undistorting "in place" into the same cv::Mat used to view
			//stereo_left_points[i]/stereo_right_points[i] would silently overwrite this
			//class's own stored point queues -- fragile even where it happens to be
			//numerically harmless, and it very nearly hid a real bug here: F is only valid in
			//the undistorted domain, so the residual below must read the undistorted points,
			//never the original distorted ones
			std::vector<cv::Point2f> undist_l, undist_r;
			cv::undistortPoints(stereo_left_points[i], undist_l, camera_matrix_l, dist_l, cv::Mat(), camera_matrix_l);
			cv::undistortPoints(stereo_right_points[i], undist_r, camera_matrix_r, dist_r, cv::Mat(), camera_matrix_r);

			std::vector<cv::Vec3f> lines_l, lines_r;
			cv::computeCorrespondEpilines(undist_l, 1, F, lines_l);
			cv::computeCorrespondEpilines(undist_r, 2, F, lines_r);

			double pair_error = 0.0;
			int npt = (int)stereo_left_points[i].size();
			for (int j = 0; j < npt; j++)
			{
				double err = std::fabs(undist_l[j].x * lines_r[j][0] + undist_l[j].y * lines_r[j][1] + lines_r[j][2])
					+ std::fabs(undist_r[j].x * lines_l[j][0] + undist_r[j].y * lines_l[j][1] + lines_l[j][2]);
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
