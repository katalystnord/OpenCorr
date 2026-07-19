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

#ifndef _CRACK_RESIDUAL_H_
#define _CRACK_RESIDUAL_H_

#include "oc_cubic_bspline.h"
#include "oc_nearest_neighbor.h"
#include "oc_shape.h"

namespace opencorr
{
	//Full-field self-consistency diagnostic, productized from DICe's
	//Crack_Locator_Post_Processor concept (dicengine/dice, src/core/DICe_PostProcessor.cpp,
	//BSD-3-Clause) -- DICe's own version is unfinished/prototype code (no machine-readable
	//output, hardcoded 8-bit PNG dumps, an undocumented one-sided residual clamp), so this is
	//a from-scratch implementation of the same underlying idea against OpenCorr's own
	//conventions, not a translation.
	//
	//What it does and why it's not redundant with StatusFlag/Uncertainty2D/SpeckleQualityMap:
	//those are all evaluated per-POI, in isolation -- two POIs individually straddling a real
	//crack can each report a perfectly good correlation (high ZNCC, good sigma/beta, good
	//speckle texture) with nothing to flag the discontinuity between them, since none of those
	//metrics look at agreement BETWEEN neighboring POIs. This diagnostic instead densely
	//re-evaluates a local displacement fit at every pixel (the same kind of KD-tree-radius-
	//search + local [1,dx,dy] least-squares fit Strain::compute() already does, see
	//oc_strain.cpp, but evaluated at every pixel and keeping the fit's intercept -- the
	//predicted inverse displacement at the query pixel itself -- instead of its gradient),
	//uses that fit to backward-warp the reference image, and diffs the warp against the real
	//target image. A sharp real discontinuity shows up as a residual spike along it, because
	//no smooth local fit spanning both sides of a genuine crack can predict both sides
	//correctly at once -- exactly the class of failure the pointwise diagnostics can't see.
	//
	//The KD-tree is built over POIs at their DEFORMED positions (poi.x + u, poi.y + v), not
	//their reference positions like Strain's own KD-tree -- because the query points here are
	//target-image pixel coordinates, and the fit needs to be local in the coordinate space the
	//query point actually lives in (Strain instead evaluates its fit exactly at a POI's own
	//known reference position, so it has no such requirement).
	//
	//Uses residual_map(y,x) == -1.f as the "not computed" sentinel (matching this codebase's
	//existing sigma/beta convention, oc_poi.h) for: too few neighbor POIs, a degenerate local
	//fit (checked via the fit's own QR rank, matching the same guard just added to
	//Strain::compute()), or the predicted source location falling outside the reference
	//image. A real residual is always >= 0 (it's an absolute difference), so -1 can't collide
	//with a genuine computed value.
	class CrackResidual2D
	{
	public:
		CrackResidual2D(float search_radius, int neighbor_number_min, int thread_number);
		~CrackResidual2D();

		void prepare(std::vector<POI2D>& poi_queue);

		//roi: if non-null, restrict the dense evaluation to pixels this shape contains
		//(reusing Shape2D::contains(), oc_shape.h) instead of the whole image -- avoids paying
		//for a per-pixel radius search + fit outside the region anyone actually cares about
		void compute(Image2D& ref_img, Image2D& tar_img, std::vector<POI2D>& poi_queue, Shape2D* roi = nullptr);

		const Eigen::MatrixXf& residualMap() const;

	private:
		float search_radius;
		int neighbor_number_min;
		int thread_number;

		Eigen::MatrixXf residual_map;

		std::vector<NearestNeighbor*> instance_pool;
		NearestNeighbor* getInstance(int tid);
	};

}//namespace opencorr

#endif //_CRACK_RESIDUAL_H_
