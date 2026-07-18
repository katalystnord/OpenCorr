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

}//namespace opencorr
