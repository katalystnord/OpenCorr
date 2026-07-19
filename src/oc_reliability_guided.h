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
 * Portions of this file are ported from ncorr_2D_cpp
 * (justinblaber/ncorr_2D_cpp), BSD-3-Clause. See THIRD-PARTY-LICENSES.md
 * for the required copyright notice, list of conditions, and disclaimer.
 */

#pragma once

#ifndef _RELIABILITY_GUIDED_H_
#define _RELIABILITY_GUIDED_H_

#include <vector>

#include "oc_dic.h"

namespace opencorr
{
	//Reliability-guided flood-fill initial-guess propagation, ported from
	//ncorr_2D_cpp's RGDIC()/worker_RGDIC() (justinblaber/ncorr_2D_cpp,
	//src/ncorr.cpp, BSD-3-Clause -- itself a C++ port of the original MATLAB
	//Ncorr's RG-DIC algorithm). Starting from one or more already-solved,
	//high-confidence "seed" POIs, propagates outward across a POI grid: the
	//best-matching solved point is always processed next (a priority queue
	//ordered by ZNCC, matching ncorr's own use of a priority queue ordered by
	//its correlation coefficient), each of its 4 grid neighbors gets an
	//initial guess by first-order Taylor extrapolation from the solved
	//point's own displacement AND its local displacement gradient (ux, uy,
	//vx, vy -- exactly the same quantities ncorr calls dv_dp1/dv_dp2/du_dp1/
	//du_dp2), the neighbor is refined by an existing local solver (ICGN2D1,
	//etc. -- ported via composition, not a new correlation kernel), and only
	//added to the queue itself if it clears both a ZNCC quality threshold and
	//a spatial jump-tolerance check (ncorr's cutoff_delta_disp) against its
	//own extrapolated guess.
	//
	//This is an efficient alternative to running FFTCC2D/FeatureAffine2D
	//independently for every POI: once a region has any solved neighbor, the
	//neighbor's own gradient already predicts a good starting point for the
	//local optimizer, so no fresh coarse search is needed there. It
	//complements those mechanisms rather than replacing them -- they are
	//still how the seed POI(s) themselves get their first solution.
	//
	//Scope note: ncorr's own implementation partitions the region across
	//threads with independent, redundantly-seeded flood-fills per partition
	//for both parallelism and robustness. This port is intentionally
	//single-threaded and single-seeded-per-call (the caller supplies seed
	//indices directly, e.g. from FFTCC2D+ICGN2D1 on a coarse grid) -- the
	//genuinely new capability here is the propagation mechanism itself; the
	//solver called per-POI (ICGN2D1 etc.) already provides its own
	//parallelism, and multi-region parallel flood-fill is a reasonable
	//future extension, not required for this to be useful.

	class ReliabilityGuided2D
	{
	public:
		ReliabilityGuided2D(int poi_number_x, int poi_number_y);

		float getZnccThreshold() const;
		void setZnccThreshold(float zncc_threshold); //default 0.5f

		float getDeltaDispTolerance() const;
		void setDeltaDispTolerance(float tolerance); //default 1.f, in pixels -- ncorr's cutoff_delta_disp

		//poi_queue must be a poi_number_x * poi_number_y row-major grid (poi_queue[row *
		//poi_number_x + col]); seed_indices are the flat indices of already-solved POIs
		//(poi_queue[i].result.zncc already set to a valid, STRICTLY POSITIVE value) to start
		//propagation from -- e.g. from FFTCC2D followed by ICGN2D1/ICGN2D2/ICLM2D1/etc.
		//Solvers whose own success convention is zncc==0 rather than a real correlation
		//score (FeatureAffine2D, which does feature-based matching rather than a
		//correlation-coefficient match) are NOT usable as a seed source here: zncc==0 is
		//indistinguishable from a POI that was never computed at all (POI2D::clear()'s own
		//default), so such a seed is silently skipped rather than accepted.
		//
		//solver.compute(POI2D*) is called to refine each newly-reached POI -- pass an already
		//prepare()'d ICGN2D1/ICGN2D2/ICLM2D1/2D2/SimplexMatch2D. The solver MUST populate
		//deformation.ux/uy/vx/vy on success, not just u/v: the flood-fill's Taylor
		//extrapolation reads a solved neighbor's own gradient to seed the NEXT hop, so a
		//translation-only solver (FFTCC2D) silently propagates a stale or zero gradient
		//across the whole filled region instead of erroring -- do not pass FFTCC2D here.
		//
		//Returns the number of POIs the flood-fill successfully reached and accepted
		//(excluding the seeds).
		int compute(std::vector<POI2D>& poi_queue, const std::vector<int>& seed_indices, DIC& solver);

	private:
		int poi_number_x, poi_number_y;
		float zncc_threshold;
		float delta_disp_tolerance;
	};

}//namespace opencorr

#endif //_RELIABILITY_GUIDED_H_
