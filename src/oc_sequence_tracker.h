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

#ifndef _SEQUENCE_TRACKER_H_
#define _SEQUENCE_TRACKER_H_

#include <vector>

#include "oc_dic.h"
#include "oc_image.h"

namespace opencorr
{
	//Multi-frame sequence tracking: OpenCorr is confirmed strictly pairwise/
	//path-independent by design (one compute() call, one reference/target
	//pair) with no sequence abstraction at all. This adds two things ported
	//from two different sources, since neither alone covers both:
	//
	//1. Reference-image update, ported from ncorr_2D_cpp's DIC_analysis()
	//   (justinblaber/ncorr_2D_cpp, src/ncorr.cpp, BSD-3-Clause): correlate
	//   against a fixed reference frame for as long as quality holds, then
	//   re-anchor the reference to the current frame once it doesn't --
	//   otherwise a long sequence eventually decorrelates entirely against
	//   an increasingly different original frame. ncorr's decision rule:
	//   take a percentile of the whole POI field's correlation coefficient
	//   (not just one point -- a single lucky/unlucky POI shouldn't drive a
	//   whole-field decision) and compare it to a threshold. Ported
	//   faithfully, expressed in OpenCorr's own higher-is-better ZNCC
	//   convention rather than ncorr's lower-is-better SSD-style coefficient.
	//   When the reference updates, cumulative displacement (relative to the
	//   very first frame) is "banked" and each POI's tracked position is
	//   re-expressed in the new reference frame's own coordinates -- ncorr's
	//   own "add displacements" / ROI-warping step.
	//
	//2. Jump-tolerance, ported from DICe's disp_jump_tol_ (dicengine/dice,
	//   DICe_Schema.cpp, BSD-3-Clause): a per-POI frame-to-frame maximum
	//   displacement-magnitude check, independent of the reference-update
	//   decision above -- catches an individual POI's tracking failure
	//   (losing the target, jumping to a spurious match) even on a frame
	//   where the field as a whole is tracking fine.
	//
	//The actual correlation math is unchanged: this class only orchestrates
	//repeated calls to a caller-supplied, already-constructed DIC solver
	//(ICGN2D1 etc.) across a frame sequence -- composition, not a new
	//correlation kernel.

	class SequenceTracker2D
	{
	public:
		SequenceTracker2D();

		float getJumpTolerance() const;
		void setJumpTolerance(float tolerance); //pixels, default 8 (DICe's own stated default)

		bool getReferenceUpdateEnabled() const;
		void setReferenceUpdateEnabled(bool enabled); //default false, matches ncorr's NO_UPDATE default policy

		float getUpdateZnccThreshold() const;
		void setUpdateZnccThreshold(float zncc_threshold); //default 0.9

		float getUpdatePercentile() const;
		//default 0.75: the fraction of POIs (by ZNCC) that must clear getUpdateZnccThreshold()
		//for the reference to stay fixed -- e.g. 0.75 means "at least 75% of POIs must still be
		//tracking well"; if fewer do, the reference re-anchors to the current frame. Matches
		//ncorr's own prctile_corrcoef/update_corrcoef pair (there expressed the other way round,
		//as a lower-is-better SSD-style percentile against a ceiling).
		void setUpdatePercentile(float percentile); //fraction in [0,1]

		struct FrameStatus
		{
			bool reference_updated; //true if this frame triggered a reference re-anchor
			int jump_rejected_count; //POIs whose frame-to-frame displacement exceeded the jump tolerance this frame
		};

		//images: the full sequence, images[0] is the original reference. poi_queue: POIs
		//given in images[0]'s coordinate system on input; on return, each POI's
		//deformation.u/v holds CUMULATIVE displacement relative to images[0] (matching
		//ncorr's own Lagrangian-perspective output convention), and result.zncc holds the
		//last successfully-computed frame's quality.
		//
		//poi_queue.deformation.u/v on INPUT is treated as displacement already accumulated
		//from BEFORE images[0] (e.g. resuming a longer sequence split across multiple calls)
		//-- it is added to, not overwritten by, what this call computes. Do NOT feed the
		//output of a single-pair solve (ReliabilityGuided2D, or ICGN2D1/FFTCC2D directly) on
		//images[0]->images[1] into this same poi_queue and then include that same pair in
		//images: the first tracked frame would re-solve images[0]->images[1] again from
		//scratch and add it on top of the already-present value, roughly doubling it. Either
		//start poi_queue at zero deformation and let this call solve every frame itself, or
		//if genuinely resuming after a prior solve, re-base poi.x/y to the already-tracked
		//position first and leave deformation.u/v holding only the prior cumulative amount.
		//
		//solver: an already-constructed DIC-derived solver (e.g. ICGN2D1) -- setImages()/
		//prepare() are called on it internally for each consecutive frame pair actually
		//needed (i.e. once per frame if the reference never updates, since the reference/
		//target pair changes every frame
		//regardless).  Returns one FrameStatus per frame beyond the first.
		std::vector<FrameStatus> compute(std::vector<Image2D>& images, std::vector<POI2D>& poi_queue, DIC& solver);

	private:
		float jump_tolerance;
		bool reference_update_enabled;
		float update_zncc_threshold;
		float update_percentile;
	};

}//namespace opencorr

#endif //_SEQUENCE_TRACKER_H_
