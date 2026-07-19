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
 * Portions of this file are ported from DICe (dicengine/dice), BSD-3-Clause.
 * See THIRD-PARTY-LICENSES.md for the required copyright notice, list of
 * conditions, and disclaimer.
 */

#pragma once

#ifndef _SUBSET_H_
#define _SUBSET_H_

#include <memory>
#include <set>
#include <utility>

#include "oc_image.h"
#include "oc_point.h"

namespace opencorr
{
	class Subset2D
	{
	public:
		Point2D center;
		int radius_x, radius_y;
		int height, width;
		int size;

		Eigen::MatrixXf eg_mat;

		Subset2D(Point2D center, int radius_x, int radius_y);
		~Subset2D() = default;

		void fill(Image2D* image);
		float zeroMeanNorm();

		//Dynamic obstruction/occlusion masking, ported from DICe's
		//Subset::turn_off_obstructed_pixels()/turn_on_previously_obstructed_pixels()/
		//is_obstructed_pixel() concept (dicengine/dice, src/base/DICe_Subset.h/.cpp,
		//BSD-3-Clause) -- issue #14. Coordinates are GLOBAL image pixel coordinates, not
		//subset-local: obstruction (a grip, tab, fixture, or off-plane object crossing the
		//field of view) is a property of the scene at the current frame, shared by every
		//subset/POI that happens to overlap it, not a per-subset-local concept -- matching
		//how DICe's own obstructed-pixel set works. Stored behind a shared_ptr (default:
		//each instance gets its own independent, empty set, exactly like a plain member)
		//so that callers can opt into sharing one live mask across multiple Subset2D
		//instances via shareObstructionMaskFrom() below, matching that documented
		//scene-wide semantics instead of requiring every subset/POI overlapping an
		//obstruction to be marked individually.
		//
		//Scope note -- this is deliberately phase 1 (data model) only, mirroring
		//oc_shape.h's own phase 1/phase 2 split for conformal ROI shapes: marking/querying
		//which global pixels are currently obstructed. What's explicitly NOT included here
		//is threading obstruction-awareness through fill()/zeroMeanNorm() (so an obstructed
		//pixel is actually excluded from the subset's mean/norm) or through ICGN/ICLM's
		//Hessian-build and ZNSSD loops (so it's excluded from the correlation itself) --
		//that's a separate, larger, cross-cutting effort (the same one oc_shape.h's own
		//scope note describes for conformal shapes) once there's a concrete need pulling it
		//forward.
		void markObstructed(int x, int y);
		void clearObstructed(int x, int y);
		void clearAllObstructed();
		bool isObstructed(int x, int y) const;

		//use the same live obstruction mask as source, instead of this instance's own
		//independent one -- marking/clearing a pixel through either instance (or any other
		//instance sharing the same mask) is then visible through all of them
		void shareObstructionMaskFrom(const Subset2D& source);

	private:
		std::shared_ptr<std::set<std::pair<int, int>>> obstructed_pixels = std::make_shared<std::set<std::pair<int, int>>>();
	};

	class Subset3D
	{
	public:
		Point3D center;
		int radius_x, radius_y, radius_z;
		int dim_x, dim_y, dim_z;
		int size;

		float*** vol_mat = nullptr;

		Subset3D(Point3D center, int radius_x, int radius_y, int radius_z);
		~Subset3D();

		void fill(Image3D* image);
		float zeroMeanNorm();
	};

}//namespace opencorr

#endif //_SUBSET_H_