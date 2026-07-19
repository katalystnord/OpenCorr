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

#include <cmath>
#include <queue>
#include <vector>

#include "oc_reliability_guided.h"

namespace opencorr
{
	namespace
	{
		struct RGQueueEntry
		{
			int index;
			float zncc;
		};

		//best (highest) ZNCC processed first, matching ncorr's own priority queue ordered
		//by its correlation coefficient (there, lower-is-better SSD-style, so the comparator
		//sense is inverted relative to ncorr's source, but the effect -- best match first -- is
		//the same)
		struct RGQueueLess
		{
			bool operator()(const RGQueueEntry& a, const RGQueueEntry& b) const
			{
				return a.zncc < b.zncc;
			}
		};
	}

	ReliabilityGuided2D::ReliabilityGuided2D(int poi_number_x, int poi_number_y)
		: poi_number_x(poi_number_x), poi_number_y(poi_number_y)
	{
		zncc_threshold = 0.5f;
		delta_disp_tolerance = 1.f;
	}

	float ReliabilityGuided2D::getZnccThreshold() const { return zncc_threshold; }
	void ReliabilityGuided2D::setZnccThreshold(float zncc_threshold) { this->zncc_threshold = zncc_threshold; }
	float ReliabilityGuided2D::getDeltaDispTolerance() const { return delta_disp_tolerance; }
	void ReliabilityGuided2D::setDeltaDispTolerance(float tolerance) { delta_disp_tolerance = tolerance; }

	int ReliabilityGuided2D::compute(std::vector<POI2D>& poi_queue, const std::vector<int>& seed_indices, DIC& solver)
	{
		int total = poi_number_x * poi_number_y;
		if ((int)poi_queue.size() != total)
		{
			throw std::string("ReliabilityGuided2D::compute(): poi_queue.size() ("
				+ std::to_string(poi_queue.size()) + ") does not match poi_number_x * poi_number_y ("
				+ std::to_string(total) + ")");
		}

		std::vector<bool> visited(total, false); //popped, or already queued -- never queued twice

		std::priority_queue<RGQueueEntry, std::vector<RGQueueEntry>, RGQueueLess> queue;

		for (int seed_idx : seed_indices)
		{
			if (seed_idx < 0 || seed_idx >= total || visited[seed_idx]) continue;
			if (poi_queue[seed_idx].result.zncc <= 0.f) continue; //seed must already carry a valid solution

			visited[seed_idx] = true;
			queue.push({ seed_idx, poi_queue[seed_idx].result.zncc });
		}

		int dr[4] = { -1, 1, 0, 0 };
		int dc[4] = { 0, 0, -1, 1 };
		int accepted = 0;

		while (!queue.empty())
		{
			RGQueueEntry current = queue.top();
			queue.pop();

			int row = current.index / poi_number_x;
			int col = current.index % poi_number_x;
			POI2D& cur_poi = poi_queue[current.index];

			for (int n = 0; n < 4; n++)
			{
				int nr = row + dr[n], nc = col + dc[n];
				if (nr < 0 || nr >= poi_number_y || nc < 0 || nc >= poi_number_x) continue;

				int nidx = nr * poi_number_x + nc;
				if (visited[nidx]) continue;
				visited[nidx] = true;

				POI2D& neighbor = poi_queue[nidx];

				//first-order Taylor extrapolation from the solved neighbor's own displacement
				//and its local displacement gradient (ux, uy, vx, vy) -- the actual trick that
				//makes flood-fill propagation cheaper than an independent coarse search per POI
				float dx = neighbor.x - cur_poi.x;
				float dy = neighbor.y - cur_poi.y;
				float guess_u = cur_poi.deformation.u + dx * cur_poi.deformation.ux + dy * cur_poi.deformation.uy;
				float guess_v = cur_poi.deformation.v + dx * cur_poi.deformation.vx + dy * cur_poi.deformation.vy;

				//solve on a local copy, not poi_queue[nidx] directly: if solver.compute()
				//throws (e.g. ICGN2D1's own "CPU thread ID over limit" precondition throw),
				//the caller-owned poi_queue[nidx] must be left exactly as it was found, not
				//holding a half-written guess that's indistinguishable from an unreached POI
				POI2D trial(neighbor);
				trial.deformation.u = guess_u;
				trial.deformation.v = guess_v;
				trial.deformation.ux = cur_poi.deformation.ux;
				trial.deformation.uy = cur_poi.deformation.uy;
				trial.deformation.vx = cur_poi.deformation.vx;
				trial.deformation.vy = cur_poi.deformation.vy;
				trial.result.zncc = 0.f; //clear any stale negative sentinel so the solver's own precondition check doesn't reject this POI outright

				solver.compute(&trial);

				float jump = std::sqrt((trial.deformation.u - guess_u) * (trial.deformation.u - guess_u)
					+ (trial.deformation.v - guess_v) * (trial.deformation.v - guess_v));

				neighbor = trial;

				if (neighbor.result.zncc > zncc_threshold && jump < delta_disp_tolerance)
				{
					queue.push({ nidx, neighbor.result.zncc });
					accepted++;
				}
				else if (!isFailureStatus(neighbor.result.zncc))
				{
					//only stamp STATUS_RELIABILITY_GUIDED_REJECTED when the solver DIDN'T already
					//report one of its own named failure codes (isFailureStatus(), oc_dic.h) --
					//this is what actually distinguishes "the correlation itself failed" from "it
					//succeeded (however poorly) but propagation rejected it," which call for
					//different tuning responses
					neighbor.result.zncc = (float)STATUS_RELIABILITY_GUIDED_REJECTED;
				}
			}
		}

		return accepted;
	}

}//namespace opencorr
