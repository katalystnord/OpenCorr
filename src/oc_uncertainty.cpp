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

#include "oc_uncertainty.h"

namespace opencorr
{
	std::unique_ptr<Uncertainty2D_> Uncertainty2D_::allocate(int subset_radius_x, int subset_radius_y)
	{
		Point2D subset_center(0, 0);

		std::unique_ptr<Uncertainty2D_> instance = std::make_unique<Uncertainty2D_>();
		instance->ref_subset = std::make_unique<Subset2D>(subset_center, subset_radius_x, subset_radius_y);
		instance->tar_subset = std::make_unique<Subset2D>(subset_center, subset_radius_x, subset_radius_y);

		return instance;
	}

	void Uncertainty2D_::release(std::unique_ptr<Uncertainty2D_>& instance)
	{
		instance->ref_subset.reset();
		instance->tar_subset.reset();
	}

	std::unique_ptr<Uncertainty2D_>& Uncertainty2D::getInstance(int tid)
	{
		if (tid >= (int)instance_pool.size())
		{
			throw std::string("CPU thread ID over limit");
		}

		return instance_pool[tid];
	}

	Uncertainty2D::Uncertainty2D(int subset_radius_x, int subset_radius_y, int thread_number)
		: ref_gradient(nullptr), tar_interp(nullptr)
	{
		this->subset_radius_x = subset_radius_x;
		this->subset_radius_y = subset_radius_y;
		this->thread_number = thread_number;

		self_adaptive = false;
		noise_std_dev = 0.f;
		perturbation = 0.1f;

		instance_pool.resize(thread_number);
#pragma omp parallel for
		for (int i = 0; i < thread_number; i++)
		{
			instance_pool[i] = Uncertainty2D_::allocate(subset_radius_x, subset_radius_y);
		}
	}

	Uncertainty2D::~Uncertainty2D()
	{
		for (auto& instance : instance_pool)
		{
			Uncertainty2D_::release(instance);
		}
		std::vector<std::unique_ptr<Uncertainty2D_>>().swap(instance_pool);
	}

	float Uncertainty2D::getPerturbation() const
	{
		return perturbation;
	}

	void Uncertainty2D::setPerturbation(float perturbation)
	{
		this->perturbation = perturbation;
	}

	float Uncertainty2D::noiseStdDev(Image2D& image)
	{
		//J. Immerkaer, "Fast Noise Variance Estimation," CVGIP (1996) 64(2): 300-302.
		//convolve with the fixed mask [[1,-2,1],[-2,4,-2],[1,-2,1]] and sum the
		//absolute response; this mask has zero response to any polynomial surface
		//up to 1st order, so the residual is (almost) pure noise.
		int width = image.width;
		int height = image.height;

		if (width < 3 || height < 3)
		{
			return 0.f;
		}

		const Eigen::MatrixXf& img = image.eg_mat;
		double sum_abs_response = 0.0;

		for (int r = 1; r < height - 1; r++)
		{
			for (int c = 1; c < width - 1; c++)
			{
				float response =
					img(r - 1, c - 1) - 2.f * img(r - 1, c) + img(r - 1, c + 1)
					- 2.f * img(r, c - 1) + 4.f * img(r, c) - 2.f * img(r, c + 1)
					+ img(r + 1, c - 1) - 2.f * img(r + 1, c) + img(r + 1, c + 1);

				sum_abs_response += std::fabs(response);
			}
		}

		//avoid M_PI: MSVC's <cmath> only defines it when _USE_MATH_DEFINES is set
		//before the include, which nothing in this build sets
		static const double pi = std::acos(-1.0);
		double scale = std::sqrt(0.5 * pi) / (6.0 * (width - 2) * (height - 2));
		return (float)(scale * sum_abs_response);
	}

	void Uncertainty2D::prepareRef()
	{
		if (ref_gradient != nullptr)
		{
			ref_gradient.reset();
		}

		ref_gradient = std::make_unique<Gradient2D4>(*ref_img);
		ref_gradient->getGradientX();
		ref_gradient->getGradientY();

		noise_std_dev = noiseStdDev(*ref_img);
	}

	void Uncertainty2D::prepareTar()
	{
		if (tar_interp != nullptr)
		{
			tar_interp.reset();
		}

		tar_interp = std::make_unique<BicubicBspline>(*tar_img);
		tar_interp->prepare();
	}

	void Uncertainty2D::prepare()
	{
		prepareRef();
		prepareTar();
	}

	float Uncertainty2D::znssd(Uncertainty2D_* instance, POI2D* poi, float ref_mean_norm, Deformation2D1& deformation, int subset_rx, int subset_ry)
	{
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
				instance->tar_subset->eg_mat(r, c) = tar_interp->compute(global_coor);
			}
		}

		float tar_mean_norm = instance->tar_subset->zeroMeanNorm();
		Eigen::MatrixXf error_img = instance->tar_subset->eg_mat * (ref_mean_norm / tar_mean_norm) - instance->ref_subset->eg_mat;

		return error_img.squaredNorm() / (ref_mean_norm * ref_mean_norm);
	}

	void Uncertainty2D::compute(POI2D* poi)
	{
		std::unique_ptr<Uncertainty2D_>& instance = getInstance(omp_get_thread_num());

		int subset_rx = subset_radius_x;
		int subset_ry = subset_radius_y;

		if (poi->result.zncc < 0
			|| poi->y - subset_ry < 0 || poi->x - subset_rx < 0
			|| poi->y + subset_ry > ref_img->height - 1 || poi->x + subset_rx > ref_img->width - 1)
		{
			//not computed -- explicit sentinels rather than leaving whatever POI2D::clear()
			//zero-initialized, or (if this POI2D is being reused across a sequence) a stale
			//value from a previously-successful compute() call on the same object
			poi->result.sigma = -1.f;
			poi->result.beta = 0.f;
			return;
		}

		int subset_width = 2 * subset_rx + 1;
		int subset_height = 2 * subset_ry + 1;

		//--- sigma: closed-form displacement-uncertainty estimate ---
		double sum_gx2 = 0.0, sum_gy2 = 0.0;
		for (int r = 0; r < subset_height; r++)
		{
			for (int c = 0; c < subset_width; c++)
			{
				int x_global = (int)poi->x + (c - subset_rx);
				int y_global = (int)poi->y + (r - subset_ry);
				float gx = ref_gradient->gradient_x(y_global, x_global);
				float gy = ref_gradient->gradient_y(y_global, x_global);
				sum_gx2 += (double)gx * gx;
				sum_gy2 += (double)gy * gy;
			}
		}

		double min_grad_energy = std::min(sum_gx2, sum_gy2);
		poi->result.sigma = min_grad_energy > 0.0
			? (float)std::sqrt(2.0 * noise_std_dev * noise_std_dev / min_grad_energy)
			: -1.f; //degenerate subset (e.g. uniform intensity): uncertainty cannot be estimated

		//--- beta: cost-function-conditioning probe around the converged solution ---
		instance->ref_subset->center = (Point2D)*poi;
		instance->ref_subset->fill(ref_img);
		float ref_mean_norm = instance->ref_subset->zeroMeanNorm();

		if (ref_mean_norm <= 0.f)
		{
			//degenerate (near-uniform-intensity) subset: same condition sigma already
			//flags via min_grad_energy above -- znssd() divides by ref_mean_norm and its
			//square, so proceeding here would silently produce NaN/Inf instead of a sentinel
			poi->result.beta = 0.f;
			return;
		}

		instance->tar_subset->center = (Point2D)*poi;

		Deformation2D1 converged(poi->deformation.u, poi->deformation.ux, poi->deformation.uy,
			poi->deformation.v, poi->deformation.vx, poi->deformation.vy);

		//the rotation probe perturbs uy/vx (a spatial gradient term), which induces a
		//displacement at the subset boundary of ~radius*rotation_step -- scaling it down by
		//the subset radius keeps that boundary displacement comparable in magnitude to the
		//flat `perturbation` the u/v probes induce everywhere, regardless of subset size, so
		//d_u/d_v/d_theta stay comparable to each other and beta stays comparable across POIs
		//solved with different subset radii
		float rotation_step = perturbation / (float)std::max(1, std::max(subset_rx, subset_ry));

		float d_u = 0.f, d_v = 0.f, d_theta = 0.f;
		{
			Deformation2D1 plus(converged), minus(converged);
			plus.setDeformation(converged.u + perturbation, converged.ux, converged.uy, converged.v, converged.vx, converged.vy);
			minus.setDeformation(converged.u - perturbation, converged.ux, converged.uy, converged.v, converged.vx, converged.vy);
			d_u = (znssd(instance.get(), poi, ref_mean_norm, plus, subset_rx, subset_ry)
				- znssd(instance.get(), poi, ref_mean_norm, minus, subset_rx, subset_ry)) / (2.f * perturbation);
		}
		{
			Deformation2D1 plus(converged), minus(converged);
			plus.setDeformation(converged.u, converged.ux, converged.uy, converged.v + perturbation, converged.vx, converged.vy);
			minus.setDeformation(converged.u, converged.ux, converged.uy, converged.v - perturbation, converged.vx, converged.vy);
			d_v = (znssd(instance.get(), poi, ref_mean_norm, plus, subset_rx, subset_ry)
				- znssd(instance.get(), poi, ref_mean_norm, minus, subset_rx, subset_ry)) / (2.f * perturbation);
		}
		{
			//infinitesimal rotation: antisymmetric perturbation of the shear terms (uy, vx),
			//holding u, v, and the symmetric (stretch) terms ux, vy fixed
			Deformation2D1 plus(converged), minus(converged);
			plus.setDeformation(converged.u, converged.ux, converged.uy - rotation_step, converged.v, converged.vx + rotation_step, converged.vy);
			minus.setDeformation(converged.u, converged.ux, converged.uy + rotation_step, converged.v, converged.vx - rotation_step, converged.vy);
			d_theta = (znssd(instance.get(), poi, ref_mean_norm, plus, subset_rx, subset_ry)
				- znssd(instance.get(), poi, ref_mean_norm, minus, subset_rx, subset_ry)) / (2.f * rotation_step);
		}

		poi->result.beta = std::sqrt(d_u * d_u + d_v * d_v + d_theta * d_theta);
	}

	void Uncertainty2D::compute(std::vector<POI2D>& poi_queue)
	{
		int queue_length = (int)poi_queue.size();
		//pin the parallel region to exactly thread_number threads: getInstance() indexes
		//instance_pool by omp_get_thread_num(), and unlike a solver that typically runs right
		//after the caller's own omp_set_num_threads() call, Uncertainty2D is meant to run as a
		//distinct, later pass (e.g. a separate GUI action) where the ambient OpenMP thread
		//count could plausibly no longer match what this instance was constructed with -- a
		//mismatch throws from inside the parallel region, which is undefined behavior (in
		//practice, std::terminate) rather than a catchable error
#pragma omp parallel for num_threads(thread_number)
		for (int i = 0; i < queue_length; i++)
		{
			compute(&poi_queue[i]);
		}
	}

}//namespace opencorr
