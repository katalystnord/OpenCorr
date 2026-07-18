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
#include <omp.h>

#include "oc_simplex.h"

namespace opencorr
{
	NelderMead::NelderMead(int max_iterations, float tolerance)
		: max_iterations(max_iterations), tolerance(tolerance)
	{}

	int NelderMead::getMaxIterations() const { return max_iterations; }
	float NelderMead::getTolerance() const { return tolerance; }
	void NelderMead::setMaxIterations(int max_iterations) { this->max_iterations = max_iterations; }
	void NelderMead::setTolerance(float tolerance) { this->tolerance = tolerance; }

	bool NelderMead::minimize(std::vector<float>& variables, const std::vector<float>& deltas,
		const std::function<float(const std::vector<float>&)>& objective_fn, int& iterations_used, float& final_cost)
	{
		//NelderMead is a generic, reusable minimizer (not limited to SimplexMatch2D's own
		//always-6-parameter usage), so an empty or mismatched input has to be rejected
		//explicitly here: an empty variables vector would otherwise divide by num_dofs=0 and
		//read cost[1] on a size-1 array a few lines below, both undefined behavior.
		if (variables.empty() || deltas.size() != variables.size())
		{
			throw std::string("NelderMead::minimize(): variables must be non-empty and deltas must be the same size as variables");
		}

		const int num_dofs = (int)variables.size();
		const int mpts = num_dofs + 1;
		const float tiny = 1e-10f;

		//initial simplex: the starting guess, plus one vertex per dimension offset by that
		//dimension's delta
		std::vector<std::vector<float>> points(mpts, variables);
		for (int i = 1; i < mpts; i++)
		{
			points[i][i - 1] += deltas[i - 1];
		}

		std::vector<float> cost(mpts);
		for (int i = 0; i < mpts; i++)
		{
			cost[i] = objective_fn(points[i]);
		}

		std::vector<float> column_sums(num_dofs, 0.f);
		for (int j = 0; j < num_dofs; j++)
		{
			for (int i = 0; i < mpts; i++)
			{
				column_sums[j] += points[i][j];
			}
		}

		std::vector<float> ptry(num_dofs);
		float old_rtol = 0.f;
		int iteration = 0;
		bool converged = false;

		for (iteration = 0; iteration < max_iterations; iteration++)
		{
			//find the best (ilo), worst (ihi), and second-worst (inhi) vertices
			int ilo = 0, ihi, inhi;
			ihi = cost[0] > cost[1] ? (inhi = 1, 0) : (inhi = 0, 1);
			for (int i = 0; i < mpts; i++)
			{
				if (cost[i] <= cost[ilo]) ilo = i;
				if (cost[i] > cost[ihi])
				{
					inhi = ihi;
					ihi = i;
				}
				else if (cost[i] > cost[inhi] && i != ihi)
				{
					inhi = i;
				}
			}

			float rtol = 2.f * std::fabs(cost[ihi] - cost[ilo]) / (std::fabs(cost[ihi]) + std::fabs(cost[ilo]) + tiny);
			bool stagnant = std::fabs(rtol - old_rtol) < tiny;
			old_rtol = rtol;

			if (rtol < tolerance || stagnant)
			{
				std::swap(cost[0], cost[ilo]);
				std::swap(points[0], points[ilo]);
				converged = true;
				break;
			}

			//reflect the worst point through the centroid of the rest
			float fac1 = 2.f / num_dofs, fac2 = fac1 + 1.f;
			for (int j = 0; j < num_dofs; j++)
			{
				ptry[j] = column_sums[j] * fac1 - points[ihi][j] * fac2;
			}
			float ytry = objective_fn(ptry);

			if (ytry < cost[ihi])
			{
				cost[ihi] = ytry;
				for (int j = 0; j < num_dofs; j++)
				{
					column_sums[j] += ptry[j] - points[ihi][j];
					points[ihi][j] = ptry[j];
				}
			}

			if (ytry <= cost[ilo])
			{
				//reflection gave a new best point -- try expanding further in that direction
				fac1 = -1.f / num_dofs, fac2 = fac1 - 2.f;
				for (int j = 0; j < num_dofs; j++)
				{
					ptry[j] = column_sums[j] * fac1 - points[ihi][j] * fac2;
				}
				ytry = objective_fn(ptry);

				if (ytry < cost[ihi])
				{
					cost[ihi] = ytry;
					for (int j = 0; j < num_dofs; j++)
					{
						column_sums[j] += ptry[j] - points[ihi][j];
						points[ihi][j] = ptry[j];
					}
				}
			}
			else if (ytry >= cost[inhi])
			{
				//reflection didn't help -- try contracting toward the centroid
				float ysave = cost[ihi];
				fac1 = 0.5f, fac2 = -0.5f;
				for (int j = 0; j < num_dofs; j++)
				{
					ptry[j] = column_sums[j] * fac1 - points[ihi][j] * fac2;
				}
				ytry = objective_fn(ptry);

				if (ytry < cost[ihi])
				{
					cost[ihi] = ytry;
					for (int j = 0; j < num_dofs; j++)
					{
						column_sums[j] += ptry[j] - points[ihi][j];
						points[ihi][j] = ptry[j];
					}
				}

				if (ytry >= ysave)
				{
					//contraction didn't help either -- shrink the whole simplex toward the best point
					for (int i = 0; i < mpts; i++)
					{
						if (i == ilo) continue;
						for (int j = 0; j < num_dofs; j++)
						{
							points[i][j] = 0.5f * (points[i][j] + points[ilo][j]);
						}
						cost[i] = objective_fn(points[i]);
					}
					for (int j = 0; j < num_dofs; j++)
					{
						float sum = 0.f;
						for (int i = 0; i < mpts; i++) sum += points[i][j];
						column_sums[j] = sum;
					}
				}
			}
		}

		iterations_used = iteration;
		//report the best vertex found, matching the fixed simplex-fix-up done on convergence above
		int ilo = 0;
		for (int i = 1; i < mpts; i++)
		{
			if (cost[i] < cost[ilo]) ilo = i;
		}
		variables = points[ilo];
		final_cost = cost[ilo];

		return converged;
	}

	std::unique_ptr<SimplexMatch2D_> SimplexMatch2D_::allocate(int subset_radius_x, int subset_radius_y)
	{
		Point2D subset_center(0, 0);

		std::unique_ptr<SimplexMatch2D_> instance = std::make_unique<SimplexMatch2D_>();
		instance->ref_subset = std::make_unique<Subset2D>(subset_center, subset_radius_x, subset_radius_y);
		instance->tar_subset = std::make_unique<Subset2D>(subset_center, subset_radius_x, subset_radius_y);

		return instance;
	}

	void SimplexMatch2D_::release(std::unique_ptr<SimplexMatch2D_>& instance)
	{
		instance->ref_subset.reset();
		instance->tar_subset.reset();
	}

	std::unique_ptr<SimplexMatch2D_>& SimplexMatch2D::getInstance(int tid)
	{
		if (tid >= (int)instance_pool.size())
		{
			throw std::string("CPU thread ID over limit");
		}

		return instance_pool[tid];
	}

	SimplexMatch2D::SimplexMatch2D(int subset_radius_x, int subset_radius_y, int thread_number)
		: tar_interp(nullptr)
	{
		this->subset_radius_x = subset_radius_x;
		this->subset_radius_y = subset_radius_y;
		this->thread_number = thread_number;

		self_adaptive = false;
		max_iterations = 200;
		tolerance = 1e-6f;
		delta_translation = 1.f;
		delta_shape = 0.01f;

		instance_pool.resize(thread_number);
#pragma omp parallel for
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool[i] = SimplexMatch2D_::allocate(subset_radius_x, subset_radius_y);
		}
	}

	SimplexMatch2D::~SimplexMatch2D()
	{
		for (auto& instance : instance_pool)
		{
			SimplexMatch2D_::release(instance);
		}
		std::vector<std::unique_ptr<SimplexMatch2D_>>().swap(instance_pool);
	}

	void SimplexMatch2D::setIteration(int max_iterations, float tolerance)
	{
		this->max_iterations = max_iterations;
		this->tolerance = tolerance;
	}

	void SimplexMatch2D::setDeltas(float delta_translation, float delta_shape)
	{
		this->delta_translation = delta_translation;
		this->delta_shape = delta_shape;
	}

	void SimplexMatch2D::prepareTar()
	{
		if (tar_interp != nullptr)
		{
			tar_interp.reset();
		}

		tar_interp = std::make_unique<BicubicBspline>(*tar_img);
		tar_interp->prepare();
	}

	void SimplexMatch2D::prepare()
	{
		prepareTar();
	}

	float SimplexMatch2D::znssd(SimplexMatch2D_* instance, float ref_mean_norm, const std::vector<float>& p, int subset_rx, int subset_ry)
	{
		Deformation2D1 deformation(p[0], p[1], p[2], p[3], p[4], p[5]);

		int subset_width = 2 * subset_rx + 1;
		int subset_height = 2 * subset_ry + 1;

		Point2D local_coor, warped_coor, global_coor;
		for (int r = 0; r < subset_height; r++)
		{
			for (int c = 0; c < subset_width; c++)
			{
				local_coor.x = c - subset_rx;
				local_coor.y = r - subset_ry;
				warped_coor = deformation.warp(local_coor);
				global_coor = instance->tar_subset->center + warped_coor;
				float sample = tar_interp->compute(global_coor);
				//BicubicBspline::compute() returns exactly -1.f (not NaN, not an exception) for
				//out-of-range coordinates -- a real pixel intensity is never negative, so this
				//is an unambiguous, cheap way to detect it without duplicating its bounds test.
				//If only PART of the warped window leaves the image (a deformation guess that
				//pushes some but not all sample points out of bounds), the resulting subset
				//isn't uniform, so the tar_mean_norm<=0.f guard below would never catch this on
				//its own -- without this check, the fabricated -1 values would silently blend
				//into the cost as if they were real intensities and bias the optimizer. Note:
				//ICGN2D1::compute() and Uncertainty2D::znssd() share this same exposure (same
				//BicubicBspline sentinel convention) and are not guarded here since fixing them
				//is outside this port's scope -- flagged separately for follow-up.
				if (sample < 0.f)
				{
					return 1e10f;
				}
				instance->tar_subset->eg_mat(r, c) = sample;
			}
		}

		float tar_mean_norm = instance->tar_subset->zeroMeanNorm();
		if (tar_mean_norm <= 0.f)
		{
			return 1e10f; //degenerate warped subset (e.g. warped outside the image): treat as a very poor match, not a NaN
		}

		Eigen::MatrixXf error_img = instance->tar_subset->eg_mat * (ref_mean_norm / tar_mean_norm) - instance->ref_subset->eg_mat;
		return error_img.squaredNorm() / (ref_mean_norm * ref_mean_norm);
	}

	void SimplexMatch2D::compute(POI2D* poi)
	{
		std::unique_ptr<SimplexMatch2D_>& instance = getInstance(omp_get_thread_num());

		int subset_rx = subset_radius_x;
		int subset_ry = subset_radius_y;

		if (poi->y - subset_ry < 0 || poi->x - subset_rx < 0
			|| poi->y + subset_ry > ref_img->height - 1 || poi->x + subset_rx > ref_img->width - 1)
		{
			poi->result.zncc = -3.f;
			return;
		}

		instance->ref_subset->center = (Point2D)*poi;
		instance->ref_subset->fill(ref_img);
		float ref_mean_norm = instance->ref_subset->zeroMeanNorm();

		if (ref_mean_norm <= 0.f)
		{
			poi->result.zncc = -2.f; //degenerate (near-uniform-intensity) reference subset
			return;
		}

		instance->tar_subset->center = (Point2D)*poi;

		std::vector<float> p = { poi->deformation.u, poi->deformation.ux, poi->deformation.uy,
			poi->deformation.v, poi->deformation.vx, poi->deformation.vy };
		std::vector<float> deltas = { delta_translation, delta_shape, delta_shape,
			delta_translation, delta_shape, delta_shape };

		NelderMead simplex(max_iterations, tolerance);
		int iterations_used = 0;
		float znssd_final = 0.f;
		simplex.minimize(p, deltas,
			[&](const std::vector<float>& variables) { return znssd(instance.get(), ref_mean_norm, variables, subset_rx, subset_ry); },
			iterations_used, znssd_final); //final_cost is the winning vertex's already-computed cost -- no redundant re-evaluation needed

		poi->deformation.u = p[0];
		poi->deformation.ux = p[1];
		poi->deformation.uy = p[2];
		poi->deformation.v = p[3];
		poi->deformation.vx = p[4];
		poi->deformation.vy = p[5];

		poi->result.zncc = 0.5f * (2.f - znssd_final);
		poi->result.iteration = (float)iterations_used;
		poi->result.convergence = znssd_final;
	}

	void SimplexMatch2D::compute(std::vector<POI2D>& poi_queue)
	{
		int queue_length = (int)poi_queue.size();
#pragma omp parallel for num_threads(thread_number)
		for (int i = 0; i < queue_length; i++)
		{
			compute(&poi_queue[i]);
		}
	}

}//namespace opencorr
