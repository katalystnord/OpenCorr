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
 * Portions of this file are ported from DICe (dicengine/dice) and from
 * ncorr_2D_cpp (justinblaber/ncorr_2D_cpp), both BSD-3-Clause. See
 * THIRD-PARTY-LICENSES.md for the required copyright notices, lists of
 * conditions, and disclaimers.
 */

#include <algorithm>
#include <cmath>

#include "oc_sequence_tracker.h"

namespace opencorr
{
	SequenceTracker2D::SequenceTracker2D()
	{
		jump_tolerance = 8.f; //DICe's own stated default (disp_jump_tol_)
		reference_update_enabled = false; //matches ncorr's NO_UPDATE default policy
		update_zncc_threshold = 0.9f;
		update_percentile = 0.75f;
	}

	float SequenceTracker2D::getJumpTolerance() const { return jump_tolerance; }
	void SequenceTracker2D::setJumpTolerance(float tolerance) { jump_tolerance = tolerance; }
	bool SequenceTracker2D::getReferenceUpdateEnabled() const { return reference_update_enabled; }
	void SequenceTracker2D::setReferenceUpdateEnabled(bool enabled) { reference_update_enabled = enabled; }
	float SequenceTracker2D::getUpdateZnccThreshold() const { return update_zncc_threshold; }
	void SequenceTracker2D::setUpdateZnccThreshold(float zncc_threshold) { update_zncc_threshold = zncc_threshold; }
	float SequenceTracker2D::getUpdatePercentile() const { return update_percentile; }
	void SequenceTracker2D::setUpdatePercentile(float percentile) { update_percentile = percentile; }

	std::vector<SequenceTracker2D::FrameStatus> SequenceTracker2D::compute(std::vector<Image2D>& images, std::vector<POI2D>& poi_queue, DIC& solver)
	{
		int n = (int)poi_queue.size();
		std::vector<FrameStatus> statuses;

		//ref_x/ref_y: each POI's current position, expressed in whatever frame is currently
		//the reference. cumulative_u/v: displacement already "banked" relative to images[0],
		//as of the last reference update. last_increment_*: the most recent within-current-
		//reference displacement/shape terms, reused both as next frame's initial guess (a
		//continuity assumption -- displacement between consecutive frames of a sequence is
		//typically close) and as the jump-tolerance baseline.
		std::vector<float> ref_x(n), ref_y(n), cumulative_u(n), cumulative_v(n);
		std::vector<float> last_increment_u(n, 0.f), last_increment_v(n, 0.f);
		std::vector<float> last_ux(n, 0.f), last_uy(n, 0.f), last_vx(n, 0.f), last_vy(n, 0.f);

		//false whenever last_increment_u/v is a fresh 0.f reset with no real prior
		//measurement behind it yet (before the very first tracked frame, and again for
		//any POI immediately after a reference update re-anchors it) -- the jump-tolerance
		//check below is a frame-to-frame CONTINUITY check (does this frame's motion look
		//discontinuous from the last one), not a bound on absolute displacement, so it has
		//nothing meaningful to compare against yet in either case and must not reject
		//whatever the true first-frame-against-this-reference displacement happens to be
		std::vector<bool> has_baseline(n, false);

		for (int i = 0; i < n; i++)
		{
			ref_x[i] = poi_queue[i].x;
			ref_y[i] = poi_queue[i].y;
			cumulative_u[i] = poi_queue[i].deformation.u;
			cumulative_v[i] = poi_queue[i].deformation.v;
		}

		int ref_idx = 0;
		int last_prepared_ref_idx = -1; //never equals a valid ref_idx, so frame 1 always fully prepares
		std::vector<POI2D> working;
		working.reserve(n);

		//if the solver exposes the reference/target split (SplittablePrepare2D, oc_dic.h),
		//a frame whose reference didn't just change only needs its target-side prepare
		//redone -- the reference-side prepare (e.g. ICGN's reference-image gradient maps)
		//depends only on ref_img, so redoing it against the SAME reference image every
		//frame is pure wasted work under the default NO_UPDATE policy (reference_update_
		//enabled=false), where ref_idx never changes after frame 1 at all. Solvers that
		//don't implement the split (checked once, not per-frame, since which solver this
		//is doesn't change during a single compute() call) fall back to the existing
		//always-correct combined prepare() every frame.
		SplittablePrepare2D* splittable = dynamic_cast<SplittablePrepare2D*>(&solver);

		for (int k = 1; k < (int)images.size(); k++)
		{
			working.clear();
			for (int i = 0; i < n; i++)
			{
				POI2D poi(Point2D(ref_x[i], ref_y[i]));
				poi.deformation.u = last_increment_u[i];
				poi.deformation.v = last_increment_v[i];
				poi.deformation.ux = last_ux[i];
				poi.deformation.uy = last_uy[i];
				poi.deformation.vx = last_vx[i];
				poi.deformation.vy = last_vy[i];
				working.push_back(poi);
			}

			solver.setImages(images[ref_idx], images[k]);
			if (splittable != nullptr && ref_idx == last_prepared_ref_idx)
			{
				splittable->prepareTar();
			}
			else
			{
				solver.prepare();
			}
			last_prepared_ref_idx = ref_idx;
			solver.compute(working);

			FrameStatus status;
			status.reference_updated = false;
			status.jump_rejected_count = 0;

			//tracked per-POI (not just via the aggregate jump_rejected_count) because the
			//reference-update step below must not bank a stale last_increment for a POI that
			//didn't actually produce a valid one this frame -- see that block for why
			std::vector<bool> succeeded_this_frame(n, false);

			for (int i = 0; i < n; i++)
			{
				if (working[i].result.zncc < 0.f)
				{
					//correlation failed outright for this POI this frame (StatusFlag,
					//oc_dic.h, defines STATUS_GOOD = 0: non-negative ZNCC is success --
					//unlike ReliabilityGuided2D's own seed-selection check, which excludes
					//zncc==0 for a different, POI2D::clear()-ambiguity reason, working[i]
					//here was unconditionally just solved by solver.compute() above, so a
					//zncc==0 result can only be a genuine FeatureAffine2D-style success,
					//never an uncomputed default) -- freeze its
					//previous cumulative value and leave last_increment untouched (so the
					//next frame's initial guess still comes from the last known-good state).
					//Still surface the failure to the caller (previously this left
					//poi_queue[i].result untouched, so a caller reading the result after
					//compute() returned would see stale data from an earlier successful frame
					//with no indication this frame's correlation actually failed).
					poi_queue[i].result = working[i].result;
					continue;
				}

				float jump = std::sqrt(std::pow(working[i].deformation.u - last_increment_u[i], 2)
					+ std::pow(working[i].deformation.v - last_increment_v[i], 2));

				if (has_baseline[i] && jump > jump_tolerance)
				{
					status.jump_rejected_count++;
					poi_queue[i].result = working[i].result;
					poi_queue[i].result.zncc = (float)STATUS_SEQUENCE_JUMP_REJECTED;
					continue; //don't accept this frame's result for this POI
				}

				last_increment_u[i] = working[i].deformation.u;
				last_increment_v[i] = working[i].deformation.v;
				last_ux[i] = working[i].deformation.ux;
				last_uy[i] = working[i].deformation.uy;
				last_vx[i] = working[i].deformation.vx;
				last_vy[i] = working[i].deformation.vy;
				has_baseline[i] = true;

				poi_queue[i].deformation.u = cumulative_u[i] + last_increment_u[i];
				poi_queue[i].deformation.v = cumulative_v[i] + last_increment_v[i];
				poi_queue[i].result = working[i].result;
				succeeded_this_frame[i] = true;
			}

			//reference-update decision: a percentile of the WHOLE field's ZNCC, not a single
			//point -- ported from ncorr's DIC_analysis() (prctile_corrcoef/update_corrcoef),
			//expressed in OpenCorr's higher-is-better ZNCC convention
			if (reference_update_enabled)
			{
				//every POI counts toward the trigger, including outright correlation
				//failures (their raw negative StatusFlag sentinel sorts below any realistic
				//update_zncc_threshold on its own) -- ncorr's own percentile is a quality
				//census of the WHOLE field, and dropping failures here would bias the
				//computed percentile upward exactly when many POIs failing is the clearest
				//signal that the reference has drifted too far and needs updating
				std::vector<float> zncc_values;
				zncc_values.reserve(n);
				for (int i = 0; i < n; i++)
					zncc_values.push_back(working[i].result.zncc);

				if (!zncc_values.empty())
				{
					std::sort(zncc_values.begin(), zncc_values.end());
					//the value below which (1 - update_percentile) fraction of POIs fall --
					//i.e. update_percentile fraction of POIs have ZNCC at or above this value
					int idx = (int)((1.f - update_percentile) * (zncc_values.size() - 1));
					idx = std::max(0, std::min((int)zncc_values.size() - 1, idx));
					float selected_zncc = zncc_values[idx];

					if (selected_zncc < update_zncc_threshold)
					{
						for (int i = 0; i < n; i++)
						{
							//only bank a POI's displacement into the new reference if it
							//actually produced a valid one THIS frame. Skipping this check
							//used to bank whatever last_increment_u/v happened to be left over
							//from an earlier successful frame for POIs that failed or were
							//jump-rejected this frame -- silently re-anchoring them to a wrong
							//position with no indication anything went wrong. A POI skipped
							//here simply keeps its last known-good ref_x/ref_y/cumulative
							//unchanged (an honest "assume no further motion since last
							//verified" fallback, not a wrong nonzero one) and its result.zncc
							//already reflects this frame's actual failure (see above), so the
							//caller can tell it needs reseeding rather than silently trusting
							//a corrupted position.
							if (!succeeded_this_frame[i]) continue;

							ref_x[i] += last_increment_u[i];
							ref_y[i] += last_increment_v[i];
							cumulative_u[i] += last_increment_u[i];
							cumulative_v[i] += last_increment_v[i];
							last_increment_u[i] = 0.f;
							last_increment_v[i] = 0.f;
							//no real prior increment against the NEW reference yet -- same
							//reasoning as the initial state before frame 1, see has_baseline's
							//own declaration above
							has_baseline[i] = false;
							//shape-gradient terms (last_ux/uy/vx/vy) are left as continuity
							//hints for the next frame, not reset -- the deformation itself
							//doesn't discontinue just because the reference frame does
						}
						ref_idx = k;
						status.reference_updated = true;
					}
				}
			}

			statuses.push_back(status);
		}

		return statuses;
	}

}//namespace opencorr
