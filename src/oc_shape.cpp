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

#include <cmath>

#include "oc_shape.h"

namespace opencorr
{
	namespace
	{
		//avoid M_PI: MSVC's <cmath> only defines it when _USE_MATH_DEFINES is set
		//before the include, which nothing in this build sets (same fix as elsewhere
		//in this fork, e.g. oc_uncertainty.cpp)
		const float pi = (float)std::acos(-1.0);
		const float two_pi = 2.f * pi;

		//signed angle from vector 1 to vector 2, in (-pi, pi] -- ported from DICe::angle_2d()
		float angle2D(float x1, float y1, float x2, float y2)
		{
			float theta1 = std::atan2(y1, x1);
			float theta2 = std::atan2(y2, x2);
			float dtheta = theta2 - theta1;
			while (dtheta > pi) dtheta -= two_pi;
			while (dtheta < -pi) dtheta += two_pi;
			return dtheta;
		}
	}

	std::vector<Point2D> Shape2D::getOwnedPixels() const
	{
		//Neither Circle2D nor Polygon2D's constructors bound how large a shape can be (only
		//that radius>0 / at least 3 vertices) -- without this check, a single malformed value
		//(e.g. a huge radius from bad user input in a GUI ROI tool) would make this loop
		//iterate a practically unbounded bounding box with no way for the caller to detect or
		//bound it in advance. The limit is generous (covers any realistic image many times
		//over) and computed in long long to avoid overflow on the multiplication itself.
		long long bbox_width = (long long)getMaxX() - getMinX() + 1;
		long long bbox_height = (long long)getMaxY() - getMinY() + 1;
		const long long max_bbox_pixels = 200000000LL; //200 megapixels
		if (bbox_width <= 0 || bbox_height <= 0 || bbox_width * bbox_height > max_bbox_pixels)
		{
			throw std::string("Shape2D::getOwnedPixels(): bounding box is invalid or unreasonably large ("
				+ std::to_string(bbox_width) + " x " + std::to_string(bbox_height) + ")");
		}

		std::vector<Point2D> pixels;
		for (int y = getMinY(); y <= getMaxY(); y++)
		{
			for (int x = getMinX(); x <= getMaxX(); x++)
			{
				if (contains(x, y))
				{
					pixels.push_back(Point2D((float)x, (float)y));
				}
			}
		}
		return pixels;
	}

	Polygon2D::Polygon2D(const std::vector<int>& vertex_x, const std::vector<int>& vertex_y)
		: vertex_x(vertex_x), vertex_y(vertex_y)
	{
		if (vertex_x.size() != vertex_y.size() || vertex_x.size() < 3)
		{
			throw std::string("Polygon2D: vertex_x/vertex_y must be the same size and have at least 3 points");
		}

		num_vertices = (int)vertex_x.size();

		//close the polygon by repeating the first vertex at the end, so edge i runs from
		//vertex i to vertex i+1 for every i in [0, num_vertices) without a special case for
		//the last edge
		this->vertex_x.push_back(vertex_x[0]);
		this->vertex_y.push_back(vertex_y[0]);

		min_x = max_x = vertex_x[0];
		min_y = max_y = vertex_y[0];
		for (int i = 1; i < num_vertices; i++)
		{
			min_x = std::min(min_x, vertex_x[i]);
			max_x = std::max(max_x, vertex_x[i]);
			min_y = std::min(min_y, vertex_y[i]);
			max_y = std::max(max_y, vertex_y[i]);
		}
	}

	bool Polygon2D::contains(int x, int y) const
	{
		//Check the boundary explicitly and first, with exact integer arithmetic (vertices and
		//query points are both integers here, so this is exact, no floating-point tie-break
		//involved). The winding-angle test below is only reliable for points strictly inside
		//or outside: a point exactly on an edge lands on atan2's tie-break for the
		//anti-parallel edge-vector pair at that edge, and different edges of the very same
		//simple polygon can tie-break oppositely -- verified concretely: on an axis-aligned
		//rectangle, a boundary point on one edge sums to ~0 (excluded) while a boundary point
		//on another edge sums to ~2*pi (included), an inconsistency with no geometric meaning.
		//Boundary points are classified as inside (a standard, deterministic convention).
		for (int i = 0; i < num_vertices; i++)
		{
			long long ex = vertex_x[i + 1] - vertex_x[i], ey = vertex_y[i + 1] - vertex_y[i];
			long long len2 = ex * ex + ey * ey;
			if (len2 == 0) continue; //degenerate (zero-length, e.g. a duplicated-consecutive
				//vertex) edge carries no boundary information -- without this, cross and dot
				//both collapse to 0 for EVERY query point, and 0>=0 && 0<=0 is unconditionally
				//true, making contains() return true for any (x, y) whatsoever
			long long px = x - vertex_x[i], py = y - vertex_y[i];
			long long cross = ex * py - ey * px;
			if (cross == 0)
			{
				//collinear with this edge's line -- inside iff also within the edge's span
				long long dot = px * ex + py * ey;
				if (dot >= 0 && dot <= len2) return true;
			}
		}

		//winding-angle test (DICe::Polygon::deactivate_pixels()): sum the signed angle
		//subtended by each edge as seen from (x, y); the point is inside iff that sum's
		//magnitude is (approximately) a full turn (pi here, since each edge contributes at
		//most one traversal and the loop covers the polygon once -- ported as DICe has it)
		float angle_sum = 0.f;
		for (int i = 0; i < num_vertices; i++)
		{
			float dx1 = (float)(vertex_x[i] - x), dy1 = (float)(vertex_y[i] - y);
			float dx2 = (float)(vertex_x[i + 1] - x), dy2 = (float)(vertex_y[i + 1] - y);
			angle_sum += angle2D(dx1, dy1, dx2, dy2);
		}
		return std::fabs(angle_sum) >= pi;
	}

	int Polygon2D::getMinX() const { return min_x; }
	int Polygon2D::getMaxX() const { return max_x; }
	int Polygon2D::getMinY() const { return min_y; }
	int Polygon2D::getMaxY() const { return max_y; }
	int Polygon2D::numVertices() const { return num_vertices; }

	Circle2D::Circle2D(int center_x, int center_y, float radius)
		: center_x(center_x), center_y(center_y), radius(radius)
	{
		if (radius <= 0.f)
		{
			throw std::string("Circle2D: radius must be positive");
		}
	}

	bool Circle2D::contains(int x, int y) const
	{
		float dx = (float)(x - center_x), dy = (float)(y - center_y);
		return dx * dx + dy * dy <= radius * radius;
	}

	int Circle2D::getMinX() const { return (int)std::floor(center_x - radius); }
	int Circle2D::getMaxX() const { return (int)std::ceil(center_x + radius); }
	int Circle2D::getMinY() const { return (int)std::floor(center_y - radius); }
	int Circle2D::getMaxY() const { return (int)std::ceil(center_y + radius); }

	RegionWithHoles2D::RegionWithHoles2D(std::unique_ptr<Shape2D> outer, std::vector<std::unique_ptr<Shape2D>> holes)
		: outer(std::move(outer)), holes(std::move(holes))
	{
		if (this->outer == nullptr)
		{
			throw std::string("RegionWithHoles2D: outer shape must not be null");
		}

		for (const auto& hole : this->holes)
		{
			if (hole == nullptr)
			{
				throw std::string("RegionWithHoles2D: hole shape must not be null");
			}
		}
	}

	bool RegionWithHoles2D::contains(int x, int y) const
	{
		if (!outer->contains(x, y))
		{
			return false;
		}

		for (const auto& hole : holes)
		{
			if (hole->contains(x, y))
			{
				return false;
			}
		}

		return true;
	}

	int RegionWithHoles2D::getMinX() const { return outer->getMinX(); }
	int RegionWithHoles2D::getMaxX() const { return outer->getMaxX(); }
	int RegionWithHoles2D::getMinY() const { return outer->getMinY(); }
	int RegionWithHoles2D::getMaxY() const { return outer->getMaxY(); }

}//namespace opencorr
