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

#include <algorithm>
#include <cmath>
#include <omp.h>

#include "oc_crack_residual.h"

namespace opencorr
{
	NearestNeighbor* CrackResidual2D::getInstance(int tid)
	{
		if (tid >= (int)instance_pool.size())
		{
			throw std::string("CPU thread ID over limit");
		}

		return instance_pool[tid];
	}

	CrackResidual2D::CrackResidual2D(float search_radius, int neighbor_number_min, int thread_number)
		: search_radius(search_radius), neighbor_number_min(neighbor_number_min), thread_number(thread_number)
	{
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool.push_back(new NearestNeighbor());
		}
	}

	CrackResidual2D::~CrackResidual2D()
	{
		for (auto& instance : instance_pool)
		{
			delete instance;
		}
		std::vector<NearestNeighbor*>().swap(instance_pool);
	}

	void CrackResidual2D::prepare(std::vector<POI2D>& poi_queue)
	{
		//KD-tree over POIs at their DEFORMED positions, not their reference positions like
		//Strain's own KD-tree (oc_strain.cpp) -- the query points in compute() below are
		//target-image pixel coordinates, so the local fit needs to be built in that same
		//coordinate space
		int queue_size = (int)poi_queue.size();
		std::vector<Point2D> deformed_positions(queue_size);
#pragma omp parallel for
		for (int i = 0; i < queue_size; i++)
		{
			deformed_positions[i].x = poi_queue[i].x + poi_queue[i].deformation.u;
			deformed_positions[i].y = poi_queue[i].y + poi_queue[i].deformation.v;
		}

#pragma omp parallel for
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool[i]->clear();
			instance_pool[i]->assignPoints(deformed_positions);
			instance_pool[i]->setSearchRadius(search_radius);
			instance_pool[i]->setSearchK(neighbor_number_min);
			instance_pool[i]->constructKdTree();
		}
	}

	void CrackResidual2D::compute(Image2D& ref_img, Image2D& tar_img, std::vector<POI2D>& poi_queue, Shape2D* roi)
	{
		residual_map = Eigen::MatrixXf::Constant(tar_img.height, tar_img.width, -1.f);

		BicubicBspline ref_interp(ref_img);
		ref_interp.prepare();

		int min_x = 0, max_x = tar_img.width - 1;
		int min_y = 0, max_y = tar_img.height - 1;
		if (roi != nullptr)
		{
			min_x = std::max(min_x, roi->getMinX());
			max_x = std::min(max_x, roi->getMaxX());
			min_y = std::max(min_y, roi->getMinY());
			max_y = std::min(max_y, roi->getMaxY());
		}

#pragma omp parallel for
		for (int y = min_y; y <= max_y; y++)
		{
			NearestNeighbor* neighbor_search = getInstance(omp_get_thread_num());

			for (int x = min_x; x <= max_x; x++)
			{
				if (roi != nullptr && !roi->contains(x, y))
				{
					continue; //outside the ROI: leave the -1.f "not computed" sentinel
				}

				Point3D query_point((float)x, (float)y, 0.f);

				//search the neighbor POIs in a subregion of given radius, falling back to KNN
				//if too few are found -- same two-step pattern as Strain::compute()
				std::vector<uint32_t> neighbor_idx;
				std::vector<nanoflann::ResultItem<uint32_t, float>> current_matches;
				int neighbor_num = neighbor_search->radiusSearch(query_point, current_matches);
				if (neighbor_num >= neighbor_number_min)
				{
					neighbor_idx.resize(neighbor_num);
					for (int i = 0; i < neighbor_num; i++)
					{
						neighbor_idx[i] = current_matches[i].first;
					}
				}
				else
				{
					std::vector<float> squared_distance;
					neighbor_num = neighbor_search->knnSearch(query_point, neighbor_idx, squared_distance);
				}

				if (neighbor_num < neighbor_number_min)
				{
					continue; //too few neighbors to fit: leave the -1.f sentinel
				}

				//local affine fit of INVERSE displacement (-u, -v) as a function of position
				//relative to the query pixel, in deformed-image coordinate space -- the fit's
				//intercept (row 0, since every coefficient row's dx/dy are relative to the
				//query pixel itself) is exactly the predicted inverse displacement AT (x, y)
				Eigen::MatrixXf coefficient_matrix = Eigen::MatrixXf::Zero(neighbor_num, 3);
				Eigen::VectorXf inv_u_vector(neighbor_num), inv_v_vector(neighbor_num);
				for (int i = 0; i < neighbor_num; i++)
				{
					POI2D& neighbor_poi = poi_queue[neighbor_idx[i]];
					float deformed_x = neighbor_poi.x + neighbor_poi.deformation.u;
					float deformed_y = neighbor_poi.y + neighbor_poi.deformation.v;

					coefficient_matrix(i, 0) = 1.f;
					coefficient_matrix(i, 1) = deformed_x - (float)x;
					coefficient_matrix(i, 2) = deformed_y - (float)y;

					inv_u_vector(i) = -neighbor_poi.deformation.u;
					inv_v_vector(i) = -neighbor_poi.deformation.v;
				}

				//degenerate-fit guard (matching the same guard just added to
				//Strain::compute(), oc_strain.cpp) -- e.g. neighbors nearly collinear in
				//deformed-image space
				Eigen::ColPivHouseholderQR<Eigen::MatrixXf> qr(coefficient_matrix);
				if (qr.rank() < 3)
				{
					continue;
				}

				Eigen::VectorXf inv_u_fit = qr.solve(inv_u_vector);
				Eigen::VectorXf inv_v_fit = qr.solve(inv_v_vector);

				float x_loc = (float)x + inv_u_fit(0);
				float y_loc = (float)y + inv_v_fit(0);

				Point2D source_location(x_loc, y_loc);
				float predicted = ref_interp.compute(source_location);
				if (predicted == -1.f)
				{
					continue; //predicted source location falls outside the reference image
				}

				residual_map(y, x) = std::fabs(tar_img.eg_mat(y, x) - predicted);
			}
		}
	}

	const Eigen::MatrixXf& CrackResidual2D::residualMap() const
	{
		return residual_map;
	}

}//namespace opencorr
