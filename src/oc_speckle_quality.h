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

#ifndef _SPECKLE_QUALITY_H_
#define _SPECKLE_QUALITY_H_

#include <memory>
#include <opencv2/opencv.hpp>

#include "oc_image.h"
#include "oc_shape.h"

namespace opencorr
{
	//Computes windowed speckle-pattern-quality maps (MIG, SSSIG) and whole-image SIFT
	//keypoint density/evenness, then (via AutoROI below) segments the single largest
	//speckled region into a Polygon2D. Addresses SurView's competitive-review finding that
	//auto-ROI/threshold-based segmentation is unclaimed by every DIC GUI reviewed, open or
	//commercial.
	//
	//MIG (mean intensity gradient) and SSSIG (sum of squared subset intensity gradients) are
	//standard DIC-literature speckle-quality metrics, computed here as windowed maps over
	//Gradient2D4's (oc_gradient.h) existing 4th-order gradient, via an integral image for
	//O(1) per-window queries rather than a naive sliding-window loop.
	//
	//The density/evenness metric is inspired by the DEF (density and evenness of features)
	//indicator in Zhou, Zuo, Jiang et al., "Image feature based quality assessment of speckle
	//patterns for digital image correlation measurement," Measurement 219 (2023) 113325,
	//https://doi.org/10.1016/j.measurement.2023.113325 -- notably co-authored by this
	//library's own maintainer. The paper's exact formula wasn't accessible while writing this
	//(SSRN/ScienceDirect both paywalled at the time), so density/evenness here is an
	//OpenCorr-original formulation of the same underlying idea (SIFT keypoint count per unit
	//area; spatial uniformity via a chi-square-style statistic over a coarse grid), not a
	//literal port -- flagged honestly rather than claimed as an exact reproduction. Uses
	//cv::SIFT::create()->detect() directly rather than OpenCorr's own SIFT2D class
	//(oc_sift.h), whose public API is shaped around matching keypoints between a reference/
	//target image pair, not scoring density on a single image.
	//
	//Scope note -- deliberately NOT implemented here: Gagnon & Day, "The Practicality of
	//Quality Assessment Metrics for Millimetre-Scale Digital Image Correlation Speckle
	//Patterns," Strain 61 (2025) e12491, https://doi.org/10.1111/str.12491, tested 10
	//published scalar speckle-quality metrics (including MIG and SSSIG) against actual
	//measured strain error and found none of them substantially agreed with real accuracy;
	//their recommendation is a synthetic-deformation self-test (warp the reference image by a
	//known field, re-correlate, measure the residual) as the more trustworthy signal. That
	//self-test loop is a materially bigger feature -- a synthetic-warp generator plus wiring
	//into an existing solver (FFTCC2D/ICGN2D1) -- and is a distinct future addition, not
	//implemented here. Any caller surfacing these metrics to a user should present them as a
	//coarse, fast indicator, not a guarantee of DIC accuracy.
	//
	//Known limitation -- AutoROI always returns a single simple polygon (the largest
	//connected speckled region). It cannot represent a specimen with two or more disjoint
	//speckled patches (a second patch is silently discarded, not flagged), or a hole/void in
	//the middle of the region. RegionWithHoles2D (oc_shape.h, issue #15) is now exactly the
	//composite-shape sibling to Polygon2D this would need for the hole case, but detect()
	//itself still only ever returns a bare Polygon2D -- lifting it to actually detect and
	//return holes (not just widen the return type) is still a future change, not done here.
	//
	//Image2D keeps two representations (cv_mat, eg_mat) that are only synchronized inside
	//load(); computeGradientMaps() re-derives eg_mat from cv_mat itself (via cv::cv2eigen)
	//rather than trusting the caller to have kept them in sync, since a caller who builds an
	//Image2D via the (width, height) constructor and writes only into cv_mat (the natural
	//pattern for camera-capture/GUI-paint workflows) would otherwise silently get MIG/SSSIG
	//computed from stale zero data.

	struct SpeckleQualityMetrics
	{
		float mean_mig = 0.f; //mean gradient magnitude over the image
		float mean_sssig = 0.f; //mean of the windowed sum-of-squared-gradient map
		float sift_density = 0.f; //SIFT keypoints per 10^4 px^2 (a readable unit, not per px^2)
		float sift_evenness = 0.f; //spatial uniformity of keypoint distribution, 0 (clustered) to 1 (uniform)
	};

	class SpeckleQualityMap
	{
	public:
		//window_radius: half-size (in px) of the square window used for the MIG/SSSIG maps;
		//must be positive (throws otherwise, matching Circle2D's validation convention in
		//oc_shape.cpp -- an unvalidated negative radius would otherwise silently corrupt the
		//windowed maps or read out of bounds)
		explicit SpeckleQualityMap(int window_radius = 15);

		//computes only the windowed MIG/SSSIG maps (no SIFT pass) -- the cheap path, for
		//callers (e.g. AutoROI) that only need the gradient-based quality map, not the
		//SIFT-derived whole-image density/evenness scalars. Also re-syncs image.eg_mat from
		//image.cv_mat as a side effect (see the class-level comment on why).
		void computeGradientMaps(Image2D& image);

		//computes everything: the gradient maps above, plus SIFT-based whole-image density/
		//evenness (a comparatively expensive full-image feature-detection pass)
		void compute(Image2D& image);

		//windowed maps, same (height, width) as the source image
		const Eigen::MatrixXf& migMap() const;
		const Eigen::MatrixXf& sssigMap() const;

		const SpeckleQualityMetrics& wholeImageMetrics() const;

	private:
		int window_radius;
		Eigen::MatrixXf mig_map, sssig_map;
		SpeckleQualityMetrics whole_image;

		void computeSiftDensityEvenness(Image2D& image);
	};

	class AutoROI
	{
	public:
		//window_radius: see SpeckleQualityMap -- same validation (must be positive)
		explicit AutoROI(int window_radius = 15);

		//returns nullptr if no valid speckled region is reliably identified: a blank/
		//untextured image, an image with no genuine two-region (speckle vs. background)
		//structure (e.g. a fully-speckled frame with only smooth lighting-driven variation --
		//checked via Otsu's own separability measure, not just "did thresholding run"), or a
		//detected region too small to be a meaningful fraction of the image (guards against
		//cv::normalize manufacturing full contrast out of pure sensor noise on a near-blank
		//image)
		std::unique_ptr<Polygon2D> detect(Image2D& image);

	private:
		int window_radius;
	};

}//namespace opencorr

#endif //_SPECKLE_QUALITY_H_
