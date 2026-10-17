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

#ifndef _SESSION_ALIGNMENT_H_
#define _SESSION_ALIGNMENT_H_

#include <vector>

#include "oc_image.h"
#include "oc_poi.h"
#include "oc_point.h"
#include "oc_shape.h"

namespace opencorr
{
	//A reference and a target photographed in separate sessions: the camera has been
	//moved, refocused or replaced between them, so the target sees the scene through a
	//different scale, rotation and shift, and need not be the same size. Nothing is
	//resampled. The camera move is estimated as a similarity, a reference point is seeded
	//from it, measured on the target as it was recorded, and the answer carried back into
	//the reference frame.
	//
	//The similarity is fitted on GROUND: a region of the reference the caller states did
	//not move between sessions (the rig, a fixture, the frame). Fitted on the specimen
	//instead, it absorbs the specimen's own rigid motion and uniform dilatation, which
	//are then lost from the measurement without trace; that is the caller's to state, and
	//the reason this class takes a region rather than guessing one.
	//
	//Three steps:
	//1. SIFT2D matches, kept where the reference keypoint lies on the ground;
	//2. a RANSAC similarity through them (OpenCV's estimateAffinePartial2D), which is
	//   accurate to a few hundredths of a pixel;
	//3. a refinement: ICGN at a grid of ground subsets, each lying wholly on the ground
	//   and seeded from step 2, then a least-squares similarity through the converged
	//   centres, reweighted to reject points the similarity does not describe (a stretch
	//   of "ground" that in fact moved). This takes the fit to a few thousandths.
	//
	//⚑ SIFT2D reads Image2D::cv_mat, not eg_mat. An image built in memory with only
	//eg_mat filled gives no matches at all, silently; here that surfaces as
	//too_few_matches rather than as an identity transform.

	//x = scale * R(angle) * X + (tx, ty): a reference pixel X to where the target shows it
	struct Similarity2D
	{
		double scale = 1.0;
		double angle = 0.0; //radians, positive turning +x toward +y (clockwise on screen, y down)
		double tx = 0.0;
		double ty = 0.0;

		Point2D map(const Point2D& reference) const;
		Point2D unmap(const Point2D& target) const;

		//the starting guess for a reference POI, as if nothing but the camera had moved:
		//u, v and the first-order gradient; higher-order terms are zeroed
		void seed(POI2D& poi) const;

		//a POI measured on the target, carried back into the reference frame: u, v become
		//the material point's movement in reference pixels, and the first-order gradient
		//the specimen's own. Higher-order terms are left as measured, in target units.
		void toReference(POI2D& poi) const;
	};

	//x = H(X), the camera TILTED against the reference: a plane seen by two cameras is
	//related by a homography, and a similarity is the special case with no tilt. h is
	//row-major with h[8] = 1: x = (h0 X + h1 Y + h2) / w, y = (h3 X + h4 Y + h5) / w,
	//w = h6 X + h7 Y + 1. Only the plane is mapped exactly; what stands off it is not.
	struct Homography2D
	{
		double h[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

		static Homography2D from(const Similarity2D& similarity);

		Point2D map(const Point2D& reference) const;
		Point2D unmap(const Point2D& target) const;

		//d map / dX at X: j = { dx/dX, dx/dY, dy/dX, dy/dY }
		void jacobian(const Point2D& reference, double j[4]) const;

		//as Similarity2D's: the camera's own local stretch and turn at the point, which
		//under a tilt differ from point to point
		void seed(POI2D& poi) const;
		void toReference(POI2D& poi) const;
	};

	enum class AlignmentStatus
	{
		ok,
		no_images, //compute() called before setImages()
		too_few_matches, //fewer than min_matches SIFT matches on the ground, or RANSAC found no consensus
		too_few_ground_points, //fewer than min_refined ground subsets converged and agreed
		known_length_unmeasured, //an end of the known length did not correlate on the target
	};

	struct AlignmentResult
	{
		AlignmentStatus status = AlignmentStatus::no_images;
		Similarity2D transform; //identity unless status is ok

		//the perspective fit through the same ground subsets. Always fitted when the
		//refinement had enough to fit it, as a check on the similarity: on a camera that
		//did not tilt the two agree. When perspective was asked for (setPerspective) it is
		//THE camera move, perspective is true, and the counts, residual and per-subset
		//report below describe it rather than the similarity.
		bool perspective = false;
		bool projective_fitted = false;
		Homography2D projective;
		float projective_rms = 0.f; //px in the target, of its own kept points

		int matches = 0; //SIFT matches with the reference keypoint on the ground
		int inliers = 0; //of those, RANSAC inliers
		int refined_points = 0; //ground subsets attempted in the refinement
		int refined_used = 0; //of those, converged and kept by the final fit
		float residual_rms = 0.f; //px in the target, of the kept points about the final fit

		//every refinement subset, in the order attempted, so a caller can say WHERE the
		//ground disagreed with the fit and not only how many subsets it lost: its centre
		//in the reference, its correlation, its distance from the final fit (px in the
		//target; NaN where it did not converge) and whether the final fit kept it. A
		//subset that correlated well and was rejected far from the fit sits on ground
		//that moved. Filled whenever the refinement ran, even if it then fell short.
		std::vector<Point2D> refined_centre;
		std::vector<float> refined_zncc;
		std::vector<float> refined_residual;
		std::vector<char> refined_kept;

		//the known length, when one was given: where each end was found on the target,
		//how well it correlated, and whether the scale was taken from it
		bool known_length_used = false;
		Point2D known_length_target[2];
		float known_length_zncc[2] = { 0.f, 0.f };
	};

	class SessionAlignment2D
	{
	public:
		//subset_radius and grid_step govern the refinement's ground subsets
		SessionAlignment2D(int subset_radius, int grid_step, int thread_number);

		void setImages(Image2D& ref_img, Image2D& tar_img);
		void setGround(const Shape2D* ground); //null: the whole reference image is ground

		void setMinimums(int min_matches, int min_refined);

		//two reference points whose distance did not change between the sessions: the
		//ends of a scale bar, two fiducials. Each is correlated on the target, and the
		//camera's scale is TAKEN from the ratio of the two lengths rather than fitted, so
		//a fit on the specimen absorbs only its rigid motion and its uniform dilatation
		//is measured. Rounded to whole pixels, since a subset is read at a whole pixel.
		//An end that does not correlate refuses the alignment (known_length_unmeasured).
		void setKnownLength(const Point2D& a, const Point2D& b);
		void clearKnownLength();

		//fit the camera move as a homography rather than a similarity: for a camera tilted
		//between the sessions. Only on GROUND: fitted on the specimen it would take up every
		//uniform strain, where the similarity takes up only the dilatation. A known length
		//is not used with it, since a tilted camera has no single scale.
		void setPerspective(bool on);

		AlignmentResult compute();

	private:
		Image2D* ref_img = nullptr;
		Image2D* tar_img = nullptr;
		const Shape2D* ground = nullptr;
		int subset_radius;
		int grid_step;
		int thread_number;
		int min_matches = 8;
		int min_refined = 6;
		bool has_known_length = false;
		bool perspective = false;
		Point2D known_length[2];
	};

}//namespace opencorr

#endif //_SESSION_ALIGNMENT_H_
