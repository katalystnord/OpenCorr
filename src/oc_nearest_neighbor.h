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

#ifndef _NEAREST_NEIGHBOR_H_
#define _NEAREST_NEIGHBOR_H_

#include <memory>

#include <nanoflann.hpp>

#include "oc_poi.h"

namespace opencorr
{
	struct Point
	{
		float x, y, z;
	};

	struct PointCloud
	{
		using coord_t = float;  //the type of each coordinate

		std::vector<Point> pts;

		//return the number of points
		inline size_t kdtree_get_point_count() const
		{
			return pts.size();
		}

		//return the dim'th component of the idx'th point
		inline float kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0)
			{
				return pts[idx].x;
			}
			else if (dim == 1)
			{
				return pts[idx].y;
			}
			else
			{
				return pts[idx].z;
			}
		}

		//optional bounding-box computation: return false to default to a standard bbox computation loop
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const
		{
			return false;
		}
	};

	class NearestNeighbor
	{
	protected:
		using KdTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3 /* dim */>;

		//point_cloud and kdt_index are held together, behind a shared_ptr, so multiple
		//NearestNeighbor instances can share one already-built tree (see shareTreeFrom()
		//below) without any one of them owning/destroying data another still references.
		//kdt_index holds a reference INTO point_cloud (nanoflann's own design, not a copy),
		//so the two must always be constructed/destroyed together as this one unit.
		struct KdTreeData
		{
			PointCloud point_cloud;
			KdTree* kdt_index = nullptr;

			~KdTreeData()
			{
				delete kdt_index;
			}
		};

		std::shared_ptr<KdTreeData> tree = std::make_shared<KdTreeData>();

		float search_radius;
		int search_k;

		//per-instance query scratch, deliberately NOT shared even when the tree is:
		//nanoflann's search calls write the query coordinates in here before searching,
		//so two threads querying concurrently through the same scratch would race
		float query_coor[3] = { 0.f };

	public:
		NearestNeighbor();
		~NearestNeighbor();

		//deleted, not merely undeclared: an ordinary compiler-generated copy would
		//silently alias tree (the shared_ptr<KdTreeData>) between two instances,
		//bypassing the deliberate, visible shareTreeFrom() below -- every current
		//caller already stores this class via unique_ptr or a non-moved stack
		//local, so this closes the gap without affecting real usage
		NearestNeighbor(const NearestNeighbor&) = delete;
		NearestNeighbor& operator=(const NearestNeighbor&) = delete;

		void assignPoints(std::vector<Point2D>& point_queue);
		void assignPoints(std::vector<POI2D>& poi_queue);
		void assignPoints(std::vector<Point3D>& point_queue);
		void assignPoints(std::vector<POI3D>& poi_queue);

		float getSearchRadius() const;
		int getSearchK() const;
		void setSearchRadius(float search_radius);
		void setSearchK(int search_k);

		void constructKdTree();
		void clear(); //clear point cloud and KD_tree

		//use an already-built tree from another instance instead of building this
		//instance's own -- for read-only concurrent querying from multiple threads
		//against the same point set, where every thread's own copy of the tree would
		//otherwise be rebuilt from scratch (identical result, pure wasted work). This
		//instance's own query scratch (query_coor) stays private, so queries through it
		//remain safe to run concurrently with queries through source or any other
		//instance sharing the same tree.
		void shareTreeFrom(const NearestNeighbor& source);

		int radiusSearch(Point3D query_point, std::vector<nanoflann::ResultItem<uint32_t, float>>& matches);
		int radiusSearch(Point3D query_point, float search_radius, std::vector<nanoflann::ResultItem<uint32_t, float>>& matches);

		int knnSearch(Point3D query_point, std::vector<uint32_t>& k_neighbor_idx, std::vector<float>& k_squared_distance);
		int knnSearch(Point3D query_point, int search_k, std::vector<uint32_t>& k_neighbor_idx, std::vector<float>& k_squared_distance);
	};

}//namespace opencorr

#endif //_NEAREST_NEIGHBOR_H_