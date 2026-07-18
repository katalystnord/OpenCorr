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

#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <functional>
#include <vector>

#include "oc_cubic_bspline.h"
#include "oc_dic.h"

namespace opencorr
{
	//Generic Nelder-Mead simplex minimizer -- a derivative-free optimizer, ported
	//from DICe::Simplex::minimize() (dicengine/dice, src/core/DICe_Simplex.cpp,
	//BSD-3-Clause). Algorithm and control-flow structure are a direct port
	//(reflect/expand/contract/shrink); reworked to take the objective as a
	//std::function rather than a virtual method (DICe uses the latter to support
	//several unrelated objective types in the same class hierarchy -- OpenCorr
	//only needs one, so a callback is simpler here) and to operate on plain
	//std::vector<float> rather than Teuchos::RCP-wrapped storage.

	class NelderMead
	{
	public:
		NelderMead(int max_iterations = 1000, float tolerance = 1e-8f);

		int getMaxIterations() const;
		float getTolerance() const;
		void setMaxIterations(int max_iterations);
		void setTolerance(float tolerance);

		//minimizes objective_fn starting from variables (in/out: initial guess in,
		//converged solution out), building the initial simplex by perturbing each
		//dimension of variables by the matching entry of deltas. Returns true if
		//converged within max_iterations; iterations_used reports how many were taken.
		bool minimize(std::vector<float>& variables, const std::vector<float>& deltas,
			const std::function<float(const std::vector<float>&)>& objective_fn, int& iterations_used);

	private:
		int max_iterations;
		float tolerance;
	};

	//Gradient-free alternative to ICGN2D1, for subsets where image gradients are
	//unreliable (low contrast, small subset size) and ICGN's Gauss-Newton
	//refinement is prone to a near-singular Hessian. Slower than ICGN -- gamma
	//(ZNSSD) is evaluated many times per POI with no gradient information to guide
	//the search -- but doesn't need gradients at all, so it still converges where
	//ICGN struggles. Uses the same 6-parameter affine model as ICGN2D1
	//(Deformation2D1: u, ux, uy, v, vx, vy) so results are directly interchangeable
	//with it; note Nelder-Mead scales poorly with dimension count, so this is
	//already close to where the method's practical limit is -- a translation-only
	//or translation+rotation fallback (fewer DOF) would converge faster and more
	//reliably if the 6-parameter model turns out to be more than a given degenerate
	//subset can actually constrain.
	//
	//ZNSSD formulation (warp -> interpolate -> zero-mean-normalize -> squared
	//error / ref_mean_norm^2) matches ICGN2D1::compute() and Uncertainty2D::znssd()
	//exactly -- same cost function, different (gradient-free) optimizer.

	class SimplexMatch2D_
	{
	public:
		std::unique_ptr<Subset2D> ref_subset;
		std::unique_ptr<Subset2D> tar_subset;

		static std::unique_ptr<SimplexMatch2D_> allocate(int subset_radius_x, int subset_radius_y);
		static void release(std::unique_ptr<SimplexMatch2D_>& instance);
	};

	class SimplexMatch2D : public DIC
	{
	private:
		std::unique_ptr<Interpolation2D> tar_interp;

		int max_iterations;
		float tolerance;

		//initial-simplex step sizes: translation (u, v) in pixels, shear/stretch
		//(ux, uy, vx, vy) as small dimensionless increments -- mirrors the
		//separate "deltas" DICe's shape functions provide per parameter type
		float delta_translation;
		float delta_shape;

		std::vector<std::unique_ptr<SimplexMatch2D_>> instance_pool;
		std::unique_ptr<SimplexMatch2D_>& getInstance(int tid);

		float znssd(SimplexMatch2D_* instance, float ref_mean_norm, const std::vector<float>& p, int subset_rx, int subset_ry);

	public:
		SimplexMatch2D(int subset_radius_x, int subset_radius_y, int thread_number);
		~SimplexMatch2D();

		void setIteration(int max_iterations, float tolerance);
		void setDeltas(float delta_translation, float delta_shape);

		void prepareTar();
		void prepare();

		void compute(POI2D* poi);
		void compute(std::vector<POI2D>& poi_queue);
	};

}//namespace opencorr

#endif //_SIMPLEX_H_
