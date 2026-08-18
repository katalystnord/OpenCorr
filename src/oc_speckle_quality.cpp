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

#include "oc_speckle_quality.h"
#include "oc_gradient.h"

namespace opencorr
{
	SpeckleQualityMap::SpeckleQualityMap(int window_radius)
		: window_radius(window_radius)
	{
		if (window_radius <= 0)
		{
			throw std::string("SpeckleQualityMap: window_radius must be positive");
		}
	}

	void SpeckleQualityMap::compute(Image2D& image)
	{
		computeGradientMaps(image);
		computeSiftDensityEvenness(image);
	}

	void SpeckleQualityMap::computeGradientMaps(Image2D& image)
	{
		//Gradient2D4 (oc_gradient.h/.cpp) reads image.eg_mat exclusively; SIFT detection
		//(computeSiftDensityEvenness) reads image.cv_mat directly. Image2D only keeps these
		//two representations in sync inside load() -- re-deriving eg_mat here rather than
		//trusting the caller removes that precondition entirely (see the class-level comment
		//in oc_speckle_quality.h for why this matters)
		cv::cv2eigen(image.cv_mat, image.eg_mat);

		int height = image.height;
		int width = image.width;

		mig_map = Eigen::MatrixXf::Zero(height, width);
		sssig_map = Eigen::MatrixXf::Zero(height, width);
		whole_image.mean_mig = 0.f;
		whole_image.mean_sssig = 0.f;

		//Gradient2D4's 5-point stencil needs at least 5px in each dimension to produce any
		//valid gradient data at all (see the "valid" tracking below) -- for a smaller image
		//there is genuinely no gradient to compute, so return the (already zeroed) maps
		//rather than silently producing values indistinguishable from "no gradient found here"
		//for a different reason
		if (width < 5 || height < 5)
		{
			return;
		}

		Gradient2D4 gradient(image);
		gradient.getGradientX();
		gradient.getGradientY();

		//Gradient2D4 leaves gradient_x zero in the outer 2 columns and gradient_y zero in the
		//outer 2 rows (its stencil never reaches them) -- a pixel only has a genuine (both-
		//axis) gradient magnitude when it's at least 2px from every edge. "valid" tracks this
		//separately from the window bounds-clamping below (r0/r1/c0/c1), so the windowed MIG
		//average divides by the count of pixels with real gradient data, not just the count of
		//in-bounds window pixels -- otherwise the known-zero stencil border silently drags the
		//average down near every image edge, which matters because AutoROI's Otsu threshold
		//runs over the whole map and is sensitive to a systematic low-value perimeter band.
		Eigen::MatrixXf integral_mag = Eigen::MatrixXf::Zero(height + 1, width + 1);
		Eigen::MatrixXf integral_sq = Eigen::MatrixXf::Zero(height + 1, width + 1);
		Eigen::MatrixXf integral_valid = Eigen::MatrixXf::Zero(height + 1, width + 1);
		float valid_mag_sum = 0.f;
		float valid_count = 0.f;
		for (int r = 0; r < height; r++)
		{
			for (int c = 0; c < width; c++)
			{
				float gx = gradient.gradient_x(r, c);
				float gy = gradient.gradient_y(r, c);
				float mag = std::sqrt(gx * gx + gy * gy);
				bool is_valid = (c >= 2 && c < width - 2 && r >= 2 && r < height - 2);

				//integral images (1 row/col larger than the source: integral(r, c) = sum of
				//all source pixels with row < r and col < c) for O(1) windowed-sum queries,
				//avoiding a naive O(window_area) sum per pixel over a potentially
				//megapixel-sized map
				integral_mag(r + 1, c + 1) = mag + integral_mag(r, c + 1) + integral_mag(r + 1, c) - integral_mag(r, c);
				integral_sq(r + 1, c + 1) = mag * mag + integral_sq(r, c + 1) + integral_sq(r + 1, c) - integral_sq(r, c);
				integral_valid(r + 1, c + 1) = (is_valid ? 1.f : 0.f) + integral_valid(r, c + 1) + integral_valid(r + 1, c) - integral_valid(r, c);

				if (is_valid)
				{
					valid_mag_sum += mag;
					valid_count += 1.f;
				}
			}
		}

		auto windowSum = [&](const Eigen::MatrixXf& integral, int r0, int r1, int c0, int c1) -> float
		{
			return integral(r1 + 1, c1 + 1) - integral(r0, c1 + 1) - integral(r1 + 1, c0) + integral(r0, c0);
		};

		for (int r = 0; r < height; r++)
		{
			int r0 = std::max(0, r - window_radius), r1 = std::min(height - 1, r + window_radius);
			for (int c = 0; c < width; c++)
			{
				int c0 = std::max(0, c - window_radius), c1 = std::min(width - 1, c + window_radius);
				float window_valid_count = windowSum(integral_valid, r0, r1, c0, c1);
				//fall back to the raw (bounds-clamped) window area if no pixel in this window
				//has real gradient data -- only reachable for a small image where the whole
				//window sits inside the 2px stencil border
				float denom = window_valid_count > 0.f ? window_valid_count : (float)(r1 - r0 + 1) * (float)(c1 - c0 + 1);
				mig_map(r, c) = windowSum(integral_mag, r0, r1, c0, c1) / denom;
				sssig_map(r, c) = windowSum(integral_sq, r0, r1, c0, c1); //sum, not mean -- matches SSSIG's definition
			}
		}

		whole_image.mean_mig = valid_count > 0.f ? valid_mag_sum / valid_count : 0.f;
		whole_image.mean_sssig = sssig_map.mean();
	}

	void SpeckleQualityMap::computeSiftDensityEvenness(Image2D& image)
	{
		cv::Ptr<cv::SIFT> detector = cv::SIFT::create();
		std::vector<cv::KeyPoint> keypoints;
		detector->detect(image.cv_mat, keypoints);

		float area = (float)image.height * (float)image.width;
		whole_image.sift_density = area > 0.f ? (float)keypoints.size() * 1e4f / area : 0.f;

		if (keypoints.empty())
		{
			whole_image.sift_evenness = 0.f;
			return;
		}

		//spatial uniformity: bin keypoints into a coarse grid, compare bin occupancy to a
		//uniform distribution via a normalized chi-square-style statistic (1 = perfectly
		//uniform coverage, 0 = every keypoint in a single bin)
		const int grid_n = 8;
		std::vector<int> bin_count(grid_n * grid_n, 0);
		for (auto& kp : keypoints)
		{
			int bx = std::min(grid_n - 1, std::max(0, (int)(kp.pt.x / image.width * grid_n)));
			int by = std::min(grid_n - 1, std::max(0, (int)(kp.pt.y / image.height * grid_n)));
			bin_count[by * grid_n + bx]++;
		}

		float n = (float)keypoints.size();
		int n_bins = grid_n * grid_n;
		float expected = n / (float)n_bins;
		float chi_sq = 0.f;
		for (int count : bin_count)
		{
			float diff = (float)count - expected;
			chi_sq += diff * diff / expected;
		}
		//worst-case chi-square (all n keypoints landing in a single bin) as the normalizer,
		//so the result is bounded to [0, 1] regardless of keypoint count or grid size
		float chi_sq_max = (n - expected) * (n - expected) / expected + (float)(n_bins - 1) * expected;
		whole_image.sift_evenness = chi_sq_max > 0.f ? std::max(0.f, 1.f - chi_sq / chi_sq_max) : 1.f;
	}

	const Eigen::MatrixXf& SpeckleQualityMap::migMap() const
	{
		return mig_map;
	}

	const Eigen::MatrixXf& SpeckleQualityMap::sssigMap() const
	{
		return sssig_map;
	}

	const SpeckleQualityMetrics& SpeckleQualityMap::wholeImageMetrics() const
	{
		return whole_image;
	}

	AutoROI::AutoROI(int window_radius)
		: window_radius(window_radius)
	{
		if (window_radius <= 0)
		{
			throw std::string("AutoROI: window_radius must be positive");
		}
	}

	std::unique_ptr<Shape2D> AutoROI::detect(Image2D& image)
	{
		SpeckleQualityMap quality_map(window_radius);
		quality_map.computeGradientMaps(image); //gradient-only path -- detect() never reads the SIFT-derived metrics, so skip that pass

		int height = image.height, width = image.width;
		const Eigen::MatrixXf& mig = quality_map.migMap();

		cv::Mat mig_mat;
		cv::eigen2cv(mig, mig_mat);

		cv::Mat mig_8u;
		cv::normalize(mig_mat, mig_8u, 0, 255, cv::NORM_MINMAX, CV_8U);

		cv::Mat binary;
		double otsu_thresh = cv::threshold(mig_8u, binary, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);

		//Otsu's method always returns SOME threshold, even when the histogram has no real
		//two-class structure -- e.g. a fully-speckled image with only smooth lighting-driven
		//MIG variation, no genuine speckled-vs-background split. In that case Otsu still finds
		//a numerically-optimal but semantically-meaningless split, silently segmenting a
		//plausible-looking but wrong partial region. Compute Otsu's own separability measure
		//(between-class variance / total variance, a standard quantity in the Otsu literature,
		//bounded to [0,1]) from the same histogram, and refuse to trust a poorly-separated
		//result rather than accept whatever split was found.
		{
			std::vector<int> hist(256, 0);
			for (int r = 0; r < mig_8u.rows; r++)
			{
				for (int c = 0; c < mig_8u.cols; c++)
				{
					hist[mig_8u.at<uchar>(r, c)]++;
				}
			}

			double total = (double)mig_8u.rows * (double)mig_8u.cols;
			double sum_all = 0.0;
			for (int i = 0; i < 256; i++)
			{
				sum_all += (double)i * hist[i];
			}
			double mean_all = total > 0.0 ? sum_all / total : 0.0;

			double total_var = 0.0;
			for (int i = 0; i < 256; i++)
			{
				total_var += hist[i] * (i - mean_all) * (i - mean_all);
			}
			total_var = total > 0.0 ? total_var / total : 0.0;

			int t = (int)std::lround(otsu_thresh);
			double w0 = 0.0, sum0 = 0.0;
			for (int i = 0; i <= t && i < 256; i++)
			{
				w0 += hist[i];
				sum0 += (double)i * hist[i];
			}
			w0 = total > 0.0 ? w0 / total : 0.0;
			double w1 = 1.0 - w0;

			//Otsu separability (between-class variance over total) was the original
			//guard here, and it does not work for this purpose. It answers "how well
			//can this histogram be cut in two", which is ~0.64 for ANY unimodal
			//spread -- a single Gaussian split at its mean gives w0=w1=0.5 and class
			//means about 0.8 sigma either side -- and cv::normalize above has already
			//removed the scale, so the value is independent of the noise amplitude.
			//Against a 0.15 floor, pure sensor noise passed comfortably, and a blank
			//frame came back segmented into a plausible-looking region covering most
			//of the image. Measured, not assumed: 0.637 for unimodal noise at every
			//sigma tried.
			//
			//What the guard actually needs to ask is whether the two classes are
			//DISTINCT, not merely whether a cut exists. Class-mean separation
			//relative to the within-class spread answers that: for a unimodal
			//distribution the two halves overlap heavily and the ratio stays near 1,
			//while a genuine speckle-versus-background split puts the means many
			//within-class standard deviations apart.
			double separation = 0.0;
			if (w0 > 0.0 && w1 > 0.0)
			{
				double mean0 = sum0 / (w0 * total);
				double mean1 = (sum_all - sum0) / (w1 * total);

				double var0 = 0.0, var1 = 0.0;
				for (int i = 0; i < 256; i++)
				{
					if (hist[i] == 0) continue;
					if (i <= t) var0 += hist[i] * (i - mean0) * (i - mean0);
					else        var1 += hist[i] * (i - mean1) * (i - mean1);
				}
				var0 /= (w0 * total);
				var1 /= (w1 * total);

				double spread = std::sqrt(var0) + std::sqrt(var1);
				if (spread > 1e-6)
				{
					separation = std::abs(mean1 - mean0) / spread;
				}
			}

			//Heuristic floor, like the area fraction below -- not a tuned or validated
			//constant. Measured values it was chosen between:
			//
			//    3.07   synthetic speckled patch on a plain background (accept)
			//    1.42 - 1.48   blank frames, sigma 0.5 to 10 (refuse)
			//    1.19   examples/2d_dic/oht_cfrp_0.bmp, a real DIC coupon (refuse)
			//
			//The real photograph scoring BELOW blank noise is not an anomaly: that
			//coupon is speckled across essentially the whole frame, so there is no
			//speckle-versus-background structure for a segmentation to find, and the
			//old code duly "found" a region covering 99.5% of the image -- a bounding
			//box identical to the frame is not a segmentation, it is a refusal
			//wearing a polygon. Declining is the better answer, and the caller can
			//offer to draw the region by hand.
			//
			//Set clear of the noise figure rather than close to the speckle one: the
			//cost of refusing wrongly is a user drawing a boundary themselves, and the
			//cost of accepting wrongly is a confident region that is not there.
			const double min_separation = 2.0;
			if (separation < min_separation)
			{
				return nullptr;
			}
		}

		//morphological cleanup: open (remove small noise specks) then close (fill small gaps).
		//kernel size scales with the quality map's own window size (so the cleanup radius
		//tracks the same spatial scale the map was computed at), but is capped relative to the
		//image's own smaller dimension -- an uncapped kernel comparable to or larger than a
		//genuinely small (but valid) speckled region can erase it entirely during MORPH_OPEN's
		//erosion step, and neither that same call's paired dilation nor the following
		//MORPH_CLOSE can recover an already-empty mask
		int kernel_size = std::max(3, (window_radius / 2) | 1);
		int max_kernel = std::max(3, std::min(height, width) / 4);
		if (kernel_size > max_kernel)
		{
			kernel_size = max_kernel | 1; //keep it odd
		}
		cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(kernel_size, kernel_size));
		cv::morphologyEx(binary, binary, cv::MORPH_OPEN, kernel);
		cv::morphologyEx(binary, binary, cv::MORPH_CLOSE, kernel);

		cv::Mat labels, stats, centroids;
		int num_labels = cv::connectedComponentsWithStats(binary, labels, stats, centroids, 8, CV_32S);
		if (num_labels <= 1) //label 0 is background; no foreground component found
		{
			return nullptr;
		}

		int best_label = 1;
		int best_area = stats.at<int>(1, cv::CC_STAT_AREA);
		for (int label = 2; label < num_labels; label++)
		{
			int area = stats.at<int>(label, cv::CC_STAT_AREA);
			if (area > best_area)
			{
				best_area = area;
				best_label = label;
			}
		}

		//require the detected region to be a meaningful fraction of the image, not just "the
		//largest of possibly-tiny noise-driven blobs" -- cv::normalize(NORM_MINMAX) above
		//manufactures full [0,255] contrast out of pure sensor noise on a near-blank image
		//(it only degenerates to constant output when max==min exactly), so num_labels<=1
		//alone is not a reliable "nothing found" signal for a real (noisy) near-blank image
		double image_area = (double)height * (double)width;
		const double min_area_fraction = 0.01; //at least 1% of the image, not a tuned/validated threshold
		if ((double)best_area < min_area_fraction * image_area)
		{
			return nullptr;
		}

		cv::Mat largest_mask = (labels == best_label);

		std::vector<std::vector<cv::Point>> contours;
		cv::findContours(largest_mask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
		if (contours.empty())
		{
			return nullptr;
		}

		//a single connected-component label can still yield more than one external contour
		//in rare pixel-level topologies (e.g. a thin bridge); keep the largest, matching the
		//same "biggest region wins" policy already applied to the connected-component choice
		size_t best_contour = 0;
		double best_contour_area = cv::contourArea(contours[0]);
		for (size_t i = 1; i < contours.size(); i++)
		{
			double contour_area = cv::contourArea(contours[i]);
			if (contour_area > best_contour_area)
			{
				best_contour_area = contour_area;
				best_contour = i;
			}
		}

		std::vector<cv::Point> simplified;
		double epsilon = 0.005 * cv::arcLength(contours[best_contour], true); //~0.5% of perimeter
		cv::approxPolyDP(contours[best_contour], simplified, epsilon, true);

		if (simplified.size() < 3)
		{
			return nullptr; //degenerate contour, can't form a polygon
		}

		std::vector<int> vertex_x, vertex_y;
		vertex_x.reserve(simplified.size());
		vertex_y.reserve(simplified.size());
		for (auto& pt : simplified)
		{
			vertex_x.push_back(pt.x);
			vertex_y.push_back(pt.y);
		}

		return std::make_unique<Polygon2D>(vertex_x, vertex_y);
	}

}//namespace opencorr
