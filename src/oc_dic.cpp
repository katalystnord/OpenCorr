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

#include "oc_dic.h"

namespace opencorr
{
	//DIC
	DIC::DIC() {}

	void DIC::setImages(Image2D& ref_img, Image2D& tar_img)
	{
		this->ref_img = &ref_img;
		this->tar_img = &tar_img;
	}

	void DIC::setSubset(int radius_x, int radius_y)
	{
		subset_radius_x = radius_x;
		subset_radius_y = radius_y;
	}

	void DIC::setSelfAdaptive(bool is_self_adaptive)
	{
		self_adaptive = is_self_adaptive;
	}


	//DVC
	DVC::DVC() {}

	void DVC::setImages(Image3D& ref_img, Image3D& tar_img)
	{
		this->ref_img = &ref_img;
		this->tar_img = &tar_img;
	}

	void DVC::setSubset(int radius_x, int radius_y, int radius_z)
	{
		subset_radius_x = radius_x;
		subset_radius_y = radius_y;
		subset_radius_z = radius_z;
	}


	std::string statusDescription(float status_code)
	{
		if (status_code >= 0.f)
		{
			return "correlation succeeded";
		}

		//exact float equality against the named constants, matching
		//ReliabilityGuided2D::isSolverFailureSentinel()'s (oc_reliability_guided.cpp) own
		//convention -- these are hardcoded literal constants the solvers assign directly,
		//never a value arrived at through computation that could drift off the literal.
		//Rounding the float to a nearest-integer enum value instead (e.g. via std::lround)
		//would misclassify a real, poor-but-genuinely-computed ZNCC that happens to fall
		//near a sentinel (e.g. -0.7, from a solver that converged to a bad match) as that
		//sentinel's specific failure reason, which would simply be false for that POI.
		if (status_code == (float)STATUS_INSUFFICIENT_FEATURES)
			return "insufficient matched keypoints in subset";
		if (status_code == (float)STATUS_DEGENERATE_INPUT)
			return "degenerate input (RANSAC consensus too small, or near-uniform-intensity subset)";
		if (status_code == (float)STATUS_INVALID_SUBSET_OR_GUESS)
			return "subset out of image bounds, or invalid (NaN/out-of-range) initial guess";
		if (status_code == (float)STATUS_MAX_ITERATIONS_REACHED)
			return "iterative solver did not converge within its iteration cap";
		if (status_code == (float)STATUS_NAN_IN_RESULT)
			return "NaN in ZNCC or displacement result";
		if (status_code == (float)STATUS_RELIABILITY_GUIDED_REJECTED)
			return "rejected by reliability-guided propagation (quality or spatial jump-tolerance gate)";
		if (status_code == (float)STATUS_SEQUENCE_JUMP_REJECTED)
			return "rejected by sequence tracking (frame-to-frame jump-tolerance gate)";
		if (status_code == (float)STATUS_HESSIAN_SINGULAR)
			return "Hessian too ill-conditioned to invert reliably (ICGN only -- ICLM's Levenberg-Marquardt damping already handles this differently)";

		//a real, computed (non-sentinel) negative ZNCC -- a poor but genuine match, not a failure
		return "correlation succeeded (poor match)";
	}

	bool sortByZNCC(const POI2D& p1, const POI2D& p2)
	{
		return p1.result.zncc > p2.result.zncc;
	}

	bool sortByDistance(const KeypointIndex& kp1, const KeypointIndex& kp2)
	{
		return kp1.distance < kp2.distance;
	}

}//namespace opencorr
