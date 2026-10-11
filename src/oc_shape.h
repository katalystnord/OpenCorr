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

#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <memory>
#include <vector>

#include "oc_point.h"

namespace opencorr
{
	//Arbitrary-region ROI data model, ported from DICe::Shape/DICe::Polygon/DICe::Circle
	//(dicengine/dice, src/base/DICe_Shape.h/.cpp, BSD-3-Clause). OpenCorr today has no
	//non-rectangular region concept at all -- POIs are a bare coordinate list the caller
	//builds by hand, and Subset2D/Subset3D are always axis-aligned rectangles/cuboids.
	//
	//Scope note -- this is deliberately phase 1 (data model) only: pixel-membership
	//testing for an arbitrary polygon or circle, i.e. "does pixel (x, y) fall inside this
	//region." That alone is already directly useful for masked POI-queue construction (see
	//getOwnedPixels() below -- build a POI grid restricted to a hand-drawn ROI). A shape
	//also reaches the correlation itself, since issue #19: a ValidPixelMask2D (oc_mask.h)
	//built from one tells the ICGN solvers, FFTCC2D and Uncertainty2D which of a subset's
	//pixels to correlate over.
	//
	//Point-in-polygon test is DICe's own method (a signed winding-angle sum around the
	//polygon's edges, not the more common ray-casting test) -- ported as-is rather than
	//substituted, since it's already a correct, standard technique.

	class Shape2D
	{
	public:
		virtual ~Shape2D() = default;

		//true if pixel (x, y) falls inside this shape
		virtual bool contains(int x, int y) const = 0;

		virtual int getMinX() const = 0;
		virtual int getMaxX() const = 0;
		virtual int getMinY() const = 0;
		virtual int getMaxY() const = 0;

		//every pixel within this shape's bounding box that contains() accepts -- unlike
		//DICe, where each derived Shape re-implements this loop itself, here it's
		//implemented once against the virtual contains()/bounding-box accessors, since nothing
		//about the enumeration itself is shape-specific
		std::vector<Point2D> getOwnedPixels() const;
	};

	//A straight-sided region of an arbitrary number of vertices, listed in order (clockwise
	//or counter-clockwise) as the boundary is traversed. Coordinates are in image
	//coordinates: (0,0) is the upper-left corner, x+ right, y+ down -- matching OpenCorr's
	//own Point2D/Image2D convention.
	class Polygon2D : public Shape2D
	{
	public:
		Polygon2D(const std::vector<int>& vertex_x, const std::vector<int>& vertex_y);

		bool contains(int x, int y) const override;
		int getMinX() const override;
		int getMaxX() const override;
		int getMinY() const override;
		int getMaxY() const override;

		int numVertices() const;

		//The boundary itself, in the order it was given. Both lists are the CLOSED form --
		//the first vertex is repeated at the end, so each holds numVertices()+1 entries --
		//because that is how contains() walks the edges and copying the closed form is
		//cheaper than reconstructing it. A caller that wants the open ring should drop the
		//last entry.
		//
		//Added because a Polygon2D was otherwise write-only: a caller handed one built by
		//AutoROI (oc_speckle_quality.h) could ask whether a pixel falls inside it, and had
		//no way to read back the boundary in order to draw, export or edit it -- which is
		//exactly what a GUI consuming auto-detection needs to do with the result.
		const std::vector<int>& vertexX() const;
		const std::vector<int>& vertexY() const;

	private:
		std::vector<int> vertex_x, vertex_y; //closed: the first vertex is repeated at the end
		int num_vertices;
		int min_x, max_x, min_y, max_y;
	};

	class Circle2D : public Shape2D
	{
	public:
		Circle2D(int center_x, int center_y, float radius);

		bool contains(int x, int y) const override;
		int getMinX() const override;
		int getMaxX() const override;
		int getMinY() const override;
		int getMaxY() const override;

	private:
		int center_x, center_y;
		float radius;
	};

	//Axis-aligned ellipse, added for SurView's region tools, where an ellipse is drawn by
	//dragging out its bounding box. The centre is a float because a box an odd number of
	//pixels wide centres on a half pixel, and rounding it would shift the whole shape by
	//half a pixel toward one side. Membership is boundary-inclusive, as Circle2D's is:
	//((x - cx) / semi_x)^2 + ((y - cy) / semi_y)^2 <= 1. Not rotatable; a rotated ellipse
	//can be added when something draws one.
	class Ellipse2D : public Shape2D
	{
	public:
		Ellipse2D(float center_x, float center_y, float semi_x, float semi_y);

		bool contains(int x, int y) const override;
		int getMinX() const override;
		int getMaxX() const override;
		int getMinY() const override;
		int getMaxY() const override;

	private:
		float center_x, center_y;
		float semi_x, semi_y;
	};

	//Multiply-connected region: an outer shape with zero or more holes cut out of it,
	//composed from any Shape2D-derived shapes (a polygon with circular holes, a circle
	//with polygonal holes, etc.) -- issue tracked as "topology-aware ROI" phase 1a.
	//
	//Scope note: this is point-membership only (contains() = outer AND NOT any hole),
	//same phase-1-only spirit as this file's own scope note above -- getOwnedPixels()
	//(inherited from Shape2D) already makes this directly useful for masked POI-queue
	//construction with a hole excluded (e.g. a bolt hole or window in a specimen).
	//Deliberately NOT included, and not currently planned unless a concrete need for it
	//arises: ncorr_2D_cpp's connectivity-aware subset/strain-neighborhood clipping (its
	//`contig_subregion_generator`), which would let a subset's own interior follow the
	//region's topology right up to a hole or crack edge instead of being excluded by a
	//margin. That capability was scoped (large -- a novel connected-component/BFS
	//primitive plus threading it through Subset2D/ICGN/ICLM's hot loops, the same
	//correlation-kernel integration this file's phase 1/phase 2 split already defers)
	//and judged not worth the cost/risk against OpenCorr's already-validated solver core
	//for the value it adds beyond simply excluding a margin around the hole from the POI
	//grid -- the value is real but concentrated in near-discontinuity precision (e.g.
	//fracture-mechanics crack-tip fields), not the common case.
	class RegionWithHoles2D : public Shape2D
	{
	public:
		RegionWithHoles2D(std::unique_ptr<Shape2D> outer, std::vector<std::unique_ptr<Shape2D>> holes);

		bool contains(int x, int y) const override;
		int getMinX() const override;
		int getMaxX() const override;
		int getMinY() const override;
		int getMaxY() const override;

	private:
		std::unique_ptr<Shape2D> outer;
		std::vector<std::unique_ptr<Shape2D>> holes;
	};

	//A region built from shapes ADDED and CUT in order, where the last shape containing a
	//pixel decides, as paint does -- so a later addition can put an island back inside an
	//earlier cut, which RegionWithHoles2D cannot express. A pixel no shape contains is
	//outside. RegionWithHoles2D is the special case "add the outer, then cut each hole".
	//Added for issue #19, where a region is handed to ValidPixelMask2D as one shape.
	//
	//The bounding box is the union of the ADDED shapes' boxes: a cut can only remove.
	//Empty (no added shape) reports min > max.
	class OrderedRegion2D : public Shape2D
	{
	public:
		OrderedRegion2D() = default;

		void add(std::unique_ptr<Shape2D> shape);
		void cut(std::unique_ptr<Shape2D> shape);

		bool contains(int x, int y) const override;
		int getMinX() const override;
		int getMaxX() const override;
		int getMinY() const override;
		int getMaxY() const override;

		int numShapes() const;

	private:
		struct Step
		{
			std::unique_ptr<Shape2D> shape;
			bool cut;
		};
		std::vector<Step> steps;
	};

}//namespace opencorr

#endif //_SHAPE_H_
