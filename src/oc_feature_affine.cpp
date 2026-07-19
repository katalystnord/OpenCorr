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

#include <numeric>
#include <random>
#include <omp.h>

#include "oc_feature_affine.h"

namespace opencorr
{
	namespace
	{
		//number of affine-matrix columns (dim + 1, for the translation/homogeneous term) --
		//also doubles as the minimum consensus-set size needed to solve for the affine
		//matrix, matching FeatureAffine2D/3D's own "essential condition to solve the
		//equation" checks (3 for 2D, 4 for 3D)
		template<typename PointT> struct AffineDim;
		template<> struct AffineDim<Point2D> { static const int value = 3; };
		template<> struct AffineDim<Point3D> { static const int value = 4; };

		inline void fillAffineRow(Eigen::MatrixXf& m, int row, const Point2D& p)
		{
			m(row, 0) = p.x;
			m(row, 1) = p.y;
			m(row, 2) = 1.f;
		}

		inline void fillAffineRow(Eigen::MatrixXf& m, int row, const Point3D& p)
		{
			m(row, 0) = p.x;
			m(row, 1) = p.y;
			m(row, 2) = p.z;
			m(row, 3) = 1.f;
		}

		inline Point2D applyAffine(const Eigen::MatrixXf& affine, const Point2D& p)
		{
			return Point2D(
				p.x * affine(0, 0) + p.y * affine(1, 0) + affine(2, 0),
				p.x * affine(0, 1) + p.y * affine(1, 1) + affine(2, 1));
		}

		inline Point3D applyAffine(const Eigen::MatrixXf& affine, const Point3D& p)
		{
			return Point3D(
				p.x * affine(0, 0) + p.y * affine(1, 0) + p.z * affine(2, 0) + affine(3, 0),
				p.x * affine(0, 1) + p.y * affine(1, 1) + p.z * affine(2, 1) + affine(3, 1),
				p.x * affine(0, 2) + p.y * affine(1, 2) + p.z * affine(2, 2) + affine(3, 2));
		}

		//RANSAC-fit an affine transform mapping ref_candidates -> tar_candidates (both
		//already POI-centered local coordinates), shared by FeatureAffine2D and
		//FeatureAffine3D -- previously duplicated near-verbatim between the two (with the
		//copy-paste wy/wz bug once caught in the 3D copy). Same structure as the original:
		//repeatedly fit an affine transform from a random minimal sample, keep the largest
		//inlier consensus set seen across trials, then re-fit using that best set. Returns
		//false (degenerate) if the final consensus set is too small to solve the equation.
		template<typename PointT>
		bool ransacAffineFit(const std::vector<PointT>& ref_candidates, const std::vector<PointT>& tar_candidates,
			const RansacConfig& ransac_config, int neighbor_number_min, std::mt19937_64& gen64,
			Eigen::MatrixXf& affine_matrix, int& trial_counter, int& max_set_size)
		{
			const int dim1 = AffineDim<PointT>::value;
			int neighbor_num = (int)ref_candidates.size();

			std::vector<int> candidate_index(neighbor_num);
			std::iota(candidate_index.begin(), candidate_index.end(), 0); //ascending order to start

			trial_counter = 0;
			float location_mean_error;
			std::vector<int> max_set;

			do
			{
				//randomly select samples
				std::shuffle(candidate_index.begin(), candidate_index.end(), gen64);

				Eigen::MatrixXf ref_neighbors(ransac_config.sample_mumber, dim1);
				Eigen::MatrixXf tar_neighbors(ransac_config.sample_mumber, dim1);
				for (int j = 0; j < ransac_config.sample_mumber; j++)
				{
					fillAffineRow(ref_neighbors, j, ref_candidates[candidate_index[j]]);
					fillAffineRow(tar_neighbors, j, tar_candidates[candidate_index[j]]);
				}
				//ref * affine = tar, thus affine is the permutation of affine matrix in the paper, where Ax=x'
				Eigen::MatrixXf affine_trial = ref_neighbors.colPivHouseholderQr().solve(tar_neighbors);

				//consensus
				std::vector<int> trial_set;
				location_mean_error = 0.f;
				for (int j = 0; j < neighbor_num; j++)
				{
					PointT predicted = applyAffine(affine_trial, ref_candidates[candidate_index[j]]);
					float estimation_error = (predicted - tar_candidates[candidate_index[j]]).vectorNorm();
					//check if the error is acceptable, keep the "good" points
					if (estimation_error < ransac_config.error_threshold)
					{
						trial_set.push_back(candidate_index[j]);
						location_mean_error += estimation_error;
					}
				}
				//replace max_set with current trial_set if the latter is larger
				if (trial_set.size() > max_set.size())
				{
					max_set.assign(trial_set.begin(), trial_set.end());
				}

				trial_counter++;
				location_mean_error /= trial_set.size();
			} while (trial_counter < ransac_config.trial_number &&
				((int)max_set.size() < neighbor_number_min || location_mean_error > ransac_config.error_threshold / neighbor_number_min));

			//calculate affine matrix according to the results of consensus
			max_set_size = (int)max_set.size();
			if (max_set_size < dim1) //essential condition to solve the equation
			{
				return false;
			}

			Eigen::MatrixXf ref_neighbors(max_set_size, dim1);
			Eigen::MatrixXf tar_neighbors(max_set_size, dim1);
			for (int i = 0; i < max_set_size; i++)
			{
				fillAffineRow(ref_neighbors, i, ref_candidates[max_set[i]]);
				fillAffineRow(tar_neighbors, i, tar_candidates[max_set[i]]);
			}
			//the method of least squares
			affine_matrix = ref_neighbors.colPivHouseholderQr().solve(tar_neighbors);
			return true;
		}
	}

	//2D implementation
	std::unique_ptr<NearestNeighbor>& FeatureAffine2D::getInstance(int tid)
	{
		if (tid >= (int)instance_pool.size())
		{
			throw std::string("CPU thread ID over limit");
		}

		return instance_pool[tid];
	}

	FeatureAffine2D::FeatureAffine2D(int radius_x, int radius_y, int thread_number)
	{
		this->subset_radius_x = radius_x;
		this->subset_radius_y = radius_y;
		neighbor_search_radius = sqrt((float)(radius_x * radius_x + radius_y * radius_y));
		neighbor_number_min = 7;
		ransac_config.error_threshold = 1.5f;
		ransac_config.sample_mumber = 3;
		ransac_config.trial_number = 20;
		this->thread_number = thread_number;

		self_adaptive = false;
		subset_feature_min = 14;
		subset_radius_min = 10;

		instance_pool.resize(thread_number);
		rng_pool.reserve(thread_number);
		std::random_device rd;
		for (int i = 0; i < thread_number; i++)
		{
			rng_pool.emplace_back(rd());
		}
#pragma omp parallel for
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool[i] = std::make_unique<NearestNeighbor>();
		}
	}

	FeatureAffine2D::~FeatureAffine2D()
	{
		for (auto& instance : instance_pool)
		{
			instance.reset();
		}
		std::vector<std::unique_ptr<NearestNeighbor>>().swap(instance_pool);
	}

	RansacConfig FeatureAffine2D::getRansacConfig() const
	{
		return ransac_config;
	}

	float FeatureAffine2D::getSearchRadius() const
	{
		return neighbor_search_radius;
	}

	int FeatureAffine2D::getNeighborMin() const
	{
		return neighbor_number_min;
	}

	void FeatureAffine2D::setSearch(float neighbor_search_radius, int neighbor_number_min)
	{
		this->neighbor_search_radius = neighbor_search_radius;
		this->neighbor_number_min = neighbor_number_min;
	}

	void FeatureAffine2D::setRansacConfig(RansacConfig ransac_config)
	{
		this->ransac_config = ransac_config;
	}

	void FeatureAffine2D::setSubsetAdjustment(int feature_min, int radius_min)
	{
		subset_feature_min = feature_min;
		subset_radius_min = radius_min;
	}

	void FeatureAffine2D::setKeypointPair(std::vector<Point2D>& ref_kp, std::vector<Point2D>& tar_kp)
	{
		this->ref_kp = ref_kp;
		this->tar_kp = tar_kp;
	}

	void FeatureAffine2D::prepare()
	{
		//build the KD-tree once (instance 0), then have every other thread's instance
		//share that same read-only tree instead of independently rebuilding an identical
		//one from the same ref_kp -- constructKdTree() is O(n log n) over the same n points
		//regardless of which instance calls it, so thread_number-1 of those builds were
		//pure duplicated work. Each instance keeps its own query scratch (query_coor),
		//so concurrent querying through different instances stays thread-safe.
		instance_pool[0]->clear();
		instance_pool[0]->assignPoints(ref_kp);
		instance_pool[0]->setSearchRadius(neighbor_search_radius);
		instance_pool[0]->setSearchK(neighbor_number_min);
		instance_pool[0]->constructKdTree();

#pragma omp parallel for
		for (int i = 1; i < thread_number; i++)
		{
			instance_pool[i]->setSearchRadius(neighbor_search_radius);
			instance_pool[i]->setSearchK(neighbor_number_min);
			instance_pool[i]->shareTreeFrom(*instance_pool[0]);
		}
	}

	void FeatureAffine2D::compute(POI2D* poi)
	{
		//set instance w.r.t. thread id
		int tid = omp_get_thread_num();
		std::unique_ptr<NearestNeighbor>& neighbor_search = getInstance(tid);

		Point3D current_point(poi->x, poi->y, 0.f);
		std::vector<Point2D> ref_candidates, tar_candidates;

		int neighbor_num = 0;

		if (self_adaptive)
		{
			float x_min = ref_img->width;
			float x_max = -1.f;
			float y_min = ref_img->height;
			float y_max = -1.f;

			//search the neighbor keypoints in a region of given radius
			std::vector<uint32_t> k_neighbor_idx;
			std::vector<float> k_squared_distance;

			neighbor_num = neighbor_search->knnSearch(current_point, subset_feature_min, k_neighbor_idx, k_squared_distance);

			if (neighbor_num < ransac_config.sample_mumber)
			{
				poi->result.zncc = (float)STATUS_INSUFFICIENT_FEATURES;
				return;
			}
			else
			{
				for (int i = 0; i < neighbor_num; i++)
				{
					ref_candidates.push_back(ref_kp[k_neighbor_idx[i]]);
					tar_candidates.push_back(tar_kp[k_neighbor_idx[i]]);

					x_min = ref_candidates[i].x < x_min ? ref_candidates[i].x : x_min;
					x_max = ref_candidates[i].x > x_max ? ref_candidates[i].x : x_max;
					y_min = ref_candidates[i].y < y_min ? ref_candidates[i].y : y_min;
					y_max = ref_candidates[i].y > y_max ? ref_candidates[i].y : y_max;
				}

				//modify POI and subset size
				if (poi->x >= x_min && poi->x <= x_max && poi->y >= y_min && poi->y <= y_max)
				{
					//if the POI is within the rectangle, set subset radius according to the farthest edges of rectangle
					poi->subset_radius.x = abs(x_max - poi->x) > abs(poi->x - x_min) ? (int)abs(x_max - poi->x) : (int)abs(poi->x - x_min);
					poi->subset_radius.y = abs(y_max - poi->y) > abs(poi->y - y_min) ? (int)abs(y_max - poi->y) : (int)abs(poi->y - y_min);
				}
				else
				{
					//if the POI is out of the rectangle, set the center of rectangle as the POI
					poi->x = (int)(0.5f * (x_max + x_min));
					poi->y = (int)(0.5f * (y_max + y_min));
					poi->subset_radius.x = (int)(0.5f * (x_max - x_min));
					poi->subset_radius.y = (int)(0.5f * (y_max - y_min));
				}

				//check if the radius is below the pre-defined minimum
				poi->subset_radius.x = poi->subset_radius.x < subset_radius_min ? subset_radius_min : poi->subset_radius.x;
				poi->subset_radius.y = poi->subset_radius.y < subset_radius_min ? subset_radius_min : poi->subset_radius.y;
			}
		}
		else
		{
			//search the neighbor keypoints in a region of given radius
			std::vector<nanoflann::ResultItem<uint32_t, float>> current_matches;
			neighbor_num = neighbor_search->radiusSearch(current_point, current_matches);

			if (neighbor_num < ransac_config.sample_mumber)
			{
				poi->result.zncc = (float)STATUS_INSUFFICIENT_FEATURES;
				return;
			}
			else
			{
				ref_candidates.resize(neighbor_num);
				tar_candidates.resize(neighbor_num);

				if (neighbor_num >= neighbor_number_min)
				{
					for (int i = 0; i < neighbor_num; i++)
					{
						ref_candidates[i] = ref_kp[current_matches[i].first];
						tar_candidates[i] = tar_kp[current_matches[i].first];
					}
				}
				else //try KNN search if the obtained neighbor keypoints are not enough
				{
					std::vector<Point2D>().swap(ref_candidates);
					std::vector<Point2D>().swap(tar_candidates);

					std::vector<uint32_t> k_neighbor_idx;
					std::vector<float> k_squared_distance;

					neighbor_num = neighbor_search->knnSearch(current_point, k_neighbor_idx, k_squared_distance);

					ref_candidates.resize(neighbor_num);
					tar_candidates.resize(neighbor_num);
					for (int i = 0; i < neighbor_num; i++)
					{
						ref_candidates[i] = ref_kp[k_neighbor_idx[i]];
						tar_candidates[i] = tar_kp[k_neighbor_idx[i]];
					}
				}
			}
		}

		//convert global coordinates to the POI-centered local coordinates
		for (int i = 0; i < neighbor_num; i++)
		{
			ref_candidates[i] = ref_candidates[i] - (Point2D)*poi;
			tar_candidates[i] = tar_candidates[i] - (Point2D)*poi;
		}

		Eigen::MatrixXf affine_matrix;
		int trial_counter = 0, max_set_size = 0;
		bool solved = ransacAffineFit(ref_candidates, tar_candidates, ransac_config, neighbor_number_min,
			rng_pool[tid], affine_matrix, trial_counter, max_set_size);

		if (!solved)
		{
			poi->result.zncc = (float)STATUS_DEGENERATE_INPUT;
			return;
		}

		//calculate the 1st order deformation according to the equivalence between affine matrix and the 1st order shape function
		poi->deformation.u = affine_matrix(2, 0);
		poi->deformation.ux = affine_matrix(0, 0) - 1.f;
		poi->deformation.uy = affine_matrix(1, 0);
		poi->deformation.v = affine_matrix(2, 1);
		poi->deformation.vx = affine_matrix(0, 1);
		poi->deformation.vy = affine_matrix(1, 1) - 1.f;

		//store results of RANSAC procedure
		poi->result.iteration = (float)trial_counter;
		poi->result.feature = (float)max_set_size;

		poi->result.zncc = 0.f;
	}

	void FeatureAffine2D::compute(std::vector<POI2D>& poi_queue)
	{
		auto queue_length = poi_queue.size();
#pragma omp parallel for
		for (int i = 0; i < queue_length; i++)
		{
			compute(&poi_queue[i]);
		}
	}



	//3D implementation
	std::unique_ptr<NearestNeighbor>& FeatureAffine3D::getInstance(int tid)
	{
		if (tid >= (int)instance_pool.size())
		{
			throw std::string("CPU thread ID over limit");
		}

		return instance_pool[tid];
	}

	FeatureAffine3D::FeatureAffine3D(int radius_x, int radius_y, int radius_z, int thread_number)
	{
		this->subset_radius_x = radius_x;
		this->subset_radius_y = radius_y;
		neighbor_search_radius = sqrt((float)(radius_x * radius_x + radius_y * radius_y + radius_z * radius_z));
		neighbor_number_min = 16;
		ransac_config.error_threshold = 3.2f;
		ransac_config.sample_mumber = 4;
		ransac_config.trial_number = 32;
		this->thread_number = thread_number;

		instance_pool.resize(thread_number);
		rng_pool.reserve(thread_number);
		std::random_device rd;
		for (int i = 0; i < thread_number; i++)
		{
			rng_pool.emplace_back(rd());
		}
#pragma omp parallel for
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool[i] = std::make_unique<NearestNeighbor>();
		}
	}

	FeatureAffine3D::~FeatureAffine3D()
	{
		for (auto& instance : instance_pool)
		{
			instance.reset();
		}
		std::vector<std::unique_ptr<NearestNeighbor>>().swap(instance_pool);
	}

	RansacConfig FeatureAffine3D::getRansacConfig() const
	{
		return ransac_config;
	}

	float FeatureAffine3D::getSearchRadius() const
	{
		return neighbor_search_radius;
	}

	int FeatureAffine3D::getNeighborMin() const
	{
		return neighbor_number_min;
	}

	void FeatureAffine3D::setSearch(float neighbor_search_radius, int neighbor_number_min)
	{
		this->neighbor_search_radius = neighbor_search_radius;
		this->neighbor_number_min = neighbor_number_min;
	}

	void FeatureAffine3D::setRansacConfig(RansacConfig ransac_config)
	{
		this->ransac_config = ransac_config;
	}

	void FeatureAffine3D::setKeypointPair(std::vector<Point3D>& ref_kp, std::vector<Point3D>& tar_kp)
	{
		this->ref_kp = ref_kp;
		this->tar_kp = tar_kp;
	}

	void FeatureAffine3D::prepare()
	{
		//see FeatureAffine2D::prepare() -- same reasoning: build once, share the rest
		instance_pool[0]->clear();
		instance_pool[0]->assignPoints(ref_kp);
		instance_pool[0]->setSearchRadius(neighbor_search_radius);
		instance_pool[0]->setSearchK(neighbor_number_min);
		instance_pool[0]->constructKdTree();

#pragma omp parallel for
		for (int i = 1; i < thread_number; i++)
		{
			instance_pool[i]->setSearchRadius(neighbor_search_radius);
			instance_pool[i]->setSearchK(neighbor_number_min);
			instance_pool[i]->shareTreeFrom(*instance_pool[0]);
		}
	}

	void FeatureAffine3D::compute(POI3D* poi)
	{
		//set instance w.r.t. thread id
		int tid = omp_get_thread_num();
		std::unique_ptr<NearestNeighbor>& neighbor_search = getInstance(tid);

		Point3D current_point(poi->x, poi->y, poi->z);
		std::vector<Point3D> ref_candidates, tar_candidates;

		//search the neighbor keypoints in a region of given radius
		std::vector<nanoflann::ResultItem<uint32_t, float>> current_matches;
		int neighbor_num = neighbor_search->radiusSearch(current_point, current_matches);

		if (neighbor_num < ransac_config.sample_mumber)
		{
			poi->result.zncc = (float)STATUS_INSUFFICIENT_FEATURES;
			return;
		}

		ref_candidates.resize(neighbor_num);
		tar_candidates.resize(neighbor_num);

		if (neighbor_num >= neighbor_number_min)
		{
			for (int i = 0; i < neighbor_num; i++)
			{
				ref_candidates[i] = ref_kp[current_matches[i].first];
				tar_candidates[i] = tar_kp[current_matches[i].first];
			}
		}
		else //try KNN search if the obtained neighbor keypoints are not enough
		{
			std::vector<Point3D>().swap(ref_candidates);
			std::vector<Point3D>().swap(tar_candidates);

			std::vector<uint32_t> k_neighbor_idx;
			std::vector<float> k_squared_distance;

			neighbor_num = neighbor_search->knnSearch(current_point, k_neighbor_idx, k_squared_distance);

			ref_candidates.resize(neighbor_num);
			tar_candidates.resize(neighbor_num);
			for (int i = 0; i < neighbor_num; i++)
			{
				ref_candidates[i] = ref_kp[k_neighbor_idx[i]];
				tar_candidates[i] = tar_kp[k_neighbor_idx[i]];
			}
		}

		//convert global coordinates to the POI-centered local coordinates
		for (int i = 0; i < neighbor_num; i++)
		{
			ref_candidates[i] = ref_candidates[i] - (Point3D)*poi;
			tar_candidates[i] = tar_candidates[i] - (Point3D)*poi;
		}

		Eigen::MatrixXf affine_matrix;
		int trial_counter = 0, max_set_size = 0;
		bool solved = ransacAffineFit(ref_candidates, tar_candidates, ransac_config, neighbor_number_min,
			rng_pool[tid], affine_matrix, trial_counter, max_set_size);

		if (!solved)
		{
			poi->result.zncc = (float)STATUS_DEGENERATE_INPUT;
			return;
		}

		//calculate the 1st order deformation according to the equivalence between affine matrix and 1st order shape function
		poi->deformation.u = affine_matrix(3, 0);
		poi->deformation.ux = affine_matrix(0, 0) - 1.f;
		poi->deformation.uy = affine_matrix(1, 0);
		poi->deformation.uz = affine_matrix(2, 0);
		poi->deformation.v = affine_matrix(3, 1);
		poi->deformation.vx = affine_matrix(0, 1);
		poi->deformation.vy = affine_matrix(1, 1) - 1.f;
		poi->deformation.vz = affine_matrix(2, 1);
		poi->deformation.w = affine_matrix(3, 2);
		poi->deformation.wx = affine_matrix(0, 2);
		poi->deformation.wy = affine_matrix(1, 2);
		poi->deformation.wz = affine_matrix(2, 2) - 1.f;

		//store results of RANSAC procedure
		poi->result.iteration = (float)trial_counter;
		poi->result.feature = (float)max_set_size;

		poi->result.zncc = 0.f;
	}

	void FeatureAffine3D::compute(std::vector<POI3D>& poi_queue)
	{
		auto queue_length = poi_queue.size();
#pragma omp parallel for
		for (int i = 0; i < queue_length; i++)
		{
			compute(&poi_queue[i]);
		}
	}

}//namespace opencorr
