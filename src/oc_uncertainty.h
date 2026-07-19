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

#ifndef _UNCERTAINTY_H_
#define _UNCERTAINTY_H_

#include "oc_cubic_bspline.h"
#include "oc_dic.h"
#include "oc_gradient.h"

namespace opencorr
{
	//Post-hoc quality metrics for an already-converged POI, computed on top
	//of a solver such as ICGN2D1 rather than replacing it:
	//
	//sigma - closed-form displacement-uncertainty estimate from image noise
	//and subset gradient content (Sutton/Schreier-style DIC error
	//propagation; the same published formula independently implemented by
	//DICe's Objective::sigma(), dicengine/dice, BSD-3-Clause).
	//noiseStdDev() implements J. Immerkaer, "Fast Noise Variance
	//Estimation," CVGIP (1996) 64(2): 300-302.
	//
	//beta - cost-function-conditioning probe: reciprocal-slope sensitivity
	//of the ZNSSD cost to a small one-sided perturbation of the rigid-body
	//subspace (u, v, and an infinitesimal rotation theta, implemented as an
	//antisymmetric perturbation of uy/vx) around the converged solution,
	//holding the remaining affine shape-function terms fixed. This isolates
	//the same (u, v, theta) subspace DICe's own Objective::beta() perturbs
	//natively, adapted to OpenCorr's 6-parameter affine model, and follows
	//DICe's own formula exactly: for each axis, evaluate ZNSSD one step of
	//epsilon on either side of the minimum, take |epsilon / (gamma - gamma_0)|
	//scaled by a per-axis factor (1e-3 for u/v, 1e-1 for theta, matching
	//DICe) on each side, and average the two. The three per-axis values
	//combine via their L2 norm. Because this is a *reciprocal* slope, a
	//LARGER beta means a FLATTER cost function near the minimum, i.e. a
	//worse-conditioned, less trustworthy match -- consistent with sigma,
	//where larger also means worse: both metrics agree larger is worse.
	//
	//Neither metric requires re-running the ICGN solve: sigma only needs
	//the reference-image gradient (already computed for the Hessian) and a
	//once-per-image noise estimate; beta only needs to re-evaluate ZNSSD at
	//a handful of perturbed points around the already-converged deformation.

	class Uncertainty2D_
	{
	public:
		std::unique_ptr<Subset2D> ref_subset;
		std::unique_ptr<Subset2D> tar_subset;

		static std::unique_ptr<Uncertainty2D_> allocate(int subset_radius_x, int subset_radius_y);
		static void release(std::unique_ptr<Uncertainty2D_>& instance);
	};

	class Uncertainty2D : public DIC
	{
	private:
		std::unique_ptr<Interpolation2D> tar_interp; //interpolation for reconstructing perturbed target subsets
		std::unique_ptr<Gradient2D4> ref_gradient; //gradient of reference image, reused from the same formulation as ICGN2D1

		float noise_std_dev; //Immerkaer noise-std-dev estimate of the reference image, computed once in prepare()
		float perturbation; //finite-difference step for beta, in pixels (u, v) and radians (theta)

		std::vector<std::unique_ptr<Uncertainty2D_>> instance_pool; //pool of instances for multi-thread processing
		std::unique_ptr<Uncertainty2D_>& getInstance(int tid);

		//ZNSSD of the given POI's subset under the given (possibly perturbed) deformation;
		//ref_subset in the instance is assumed already filled and mean-normalized by the caller
		float znssd(Uncertainty2D_* instance, POI2D* poi, float ref_mean_norm, Deformation2D1& deformation, int subset_rx, int subset_ry);

	public:
		Uncertainty2D(int subset_radius_x, int subset_radius_y, int thread_number);
		~Uncertainty2D();

		float getPerturbation() const;
		void setPerturbation(float perturbation); //default 0.1, matching the DICe reference implementation this is modeled on

		void prepareRef(); //gradient map of ref image + global noise estimate
		void prepareTar(); //interpolation coefficient look-up table of tar image
		void prepare();

		void compute(POI2D* poi);
		void compute(std::vector<POI2D>& poi_queue);

		static float noiseStdDev(Image2D& image); //Immerkaer's fast noise-variance estimator
	};

}//namespace opencorr

#endif //_UNCERTAINTY_H_
