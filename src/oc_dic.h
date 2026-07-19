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

#ifndef _DIC_H_
#define _DIC_H_

#include<memory>
#include<string>

#include "oc_poi.h"
#include "oc_subset.h"

namespace opencorr
{

	//named status codes for POI2D/POI3D::result.zncc's "abnormal value" sentinel range.
	//values match the numeric codes solvers have always assigned (so most existing
	//comparisons against raw numbers keep working unchanged) with one deliberate fix:
	//NR2D1's out-of-bounds-subset/invalid-guess check previously stamped -1 (colliding
	//with FeatureAffine's unrelated "insufficient features" meaning) and used a
	//different threshold than ICGN/ICLM's own version of the identical check -- both
	//are now aligned to STATUS_INVALID_SUBSET_OR_GUESS's shared convention, see
	//oc_nr.cpp. A given code can also mean slightly different things across solver
	//classes (noted per-value below, and not perfectly consistent even among them --
	//e.g. SimplexMatch2D unconditionally overwrites on this check where ICGN/ICLM/NR
	//preserve a pre-existing negative value instead); use statusDescription() below
	//for a precise, class-agnostic human-readable string.
	enum StatusFlag : int
	{
		STATUS_GOOD = 0, //non-negative ZNCC: correlation succeeded
		STATUS_INSUFFICIENT_FEATURES = -1, //too few matched keypoints in subset (FeatureAffine)
		STATUS_DEGENERATE_INPUT = -2, //RANSAC consensus set too small (FeatureAffine), or near-uniform-intensity reference subset (SimplexMatch2D)
		STATUS_INVALID_SUBSET_OR_GUESS = -3, //subset out of image bounds, or incoming initial guess is NaN/negative/out of range (ICGN, ICLM, NR, SimplexMatch2D)
		STATUS_MAX_ITERATIONS_REACHED = -4, //iterative solver hit its iteration cap without converging (ICGN, ICLM, NR)
		STATUS_NAN_IN_RESULT = -5, //ZNCC or displacement came out NaN (ICGN, ICLM, NR)
		STATUS_RELIABILITY_GUIDED_REJECTED = -6, //quality or spatial jump-tolerance gate not met during RG-DIC propagation (ReliabilityGuided2D)
		STATUS_SEQUENCE_JUMP_REJECTED = -7, //frame-to-frame (temporal) jump-tolerance gate not met (SequenceTracker2D)
		STATUS_HESSIAN_SINGULAR = -8, //Hessian too ill-conditioned to invert reliably (ICGN only -- ICLM's Levenberg-Marquardt damping already handles this differently, see oc_iclm.cpp)
	};

	//human-readable description of a StatusFlag code (or any non-negative ZNCC, described as success).
	//takes a float since that's the storage type of Result2D::zncc/Result3D::zncc.
	std::string statusDescription(float status_code);

	//true if zncc is exactly one of the solvers' own named failure codes (StatusFlag
	//above), as opposed to a real (if possibly low or negative) computed correlation
	//value -- a solver can legitimately converge (hit none of its own failure branches)
	//yet still report a low or even negative ZNCC for a genuinely poor match, so
	//zncc>=0 is not by itself a reliable "solver succeeded" test. Exact float equality
	//is safe here since these are hardcoded literal constants the solvers assign
	//directly, never a value arrived at through computation that could drift off the
	//literal. Single source of truth shared with statusDescription() below -- both are
	//driven by the same table, so they cannot silently drift apart the way this
	//function's previous incarnation (ReliabilityGuided2D's own local
	//isSolverFailureSentinel(), which had quietly gone stale missing two of the eight
	//codes) did.
	bool isFailureStatus(float zncc);

	//structure for brute force searching
	struct KeypointIndex
	{
		int kp_idx; //index in keypoint queue
		float distance; //Euclidean distance to the POI
	};

	//Optional capability: a solver whose prepare() step is naturally separable into a
	//reference-image phase (depends only on ref_img -- expensive to redo when ref_img
	//hasn't actually changed) and a target-image phase (must be redone whenever tar_img
	//changes) can implement this alongside DIC, so a caller processing many frames
	//against a fixed reference (SequenceTracker2D) can skip the redundant reference-side
	//work on frames where the reference didn't change, instead of paying for a full
	//prepare() every single frame. Purely additive -- a solver that doesn't implement
	//this is used exactly as before, via its own combined prepare().
	class SplittablePrepare2D
	{
	public:
		virtual ~SplittablePrepare2D() = default;
		virtual void prepareRef() = 0;
		virtual void prepareTar() = 0;
	};

	class DIC
	{
	public:
		Image2D* ref_img = nullptr;
		Image2D* tar_img = nullptr;

		int subset_radius_x, subset_radius_y;
		int thread_number; //OpenMP thread number
		bool self_adaptive;

		DIC();
		virtual ~DIC() = default;

		void setImages(Image2D& ref_img, Image2D& tar_img);
		void setSubset(int radius_x, int radius_y);
		void setSelfAdaptive(bool is_self_adaptive); //select if the subset is automatically set or manually set

		virtual void prepare() = 0;
		virtual void compute(POI2D* poi) = 0;
		virtual void compute(std::vector<POI2D>& poi_queue) = 0;

	};

	class DVC
	{
	public:
		Image3D* ref_img = nullptr;
		Image3D* tar_img = nullptr;

		int subset_radius_x, subset_radius_y, subset_radius_z;
		int thread_number; //OpenMP thread number

		DVC();
		virtual ~DVC() = default;

		void setImages(Image3D& ref_img, Image3D& tar_img);
		void setSubset(int radius_x, int radius_y, int radius_z);

		virtual void prepare() = 0;
		virtual void compute(POI3D* POI) = 0;
		virtual void compute(std::vector<POI3D>& poi_queue) = 0;
	};

	bool sortByZNCC(const POI2D& p1, const POI2D& p2);

	bool sortByDistance(const KeypointIndex& kp1, const KeypointIndex& kp2);

}//namespace opencorr

#endif //_DIC_H_
