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

#ifndef _MASK_H_
#define _MASK_H_

#include <cstdint>

#include "oc_array.h"
#include "oc_shape.h"

namespace opencorr
{
	//Which pixels of the REFERENCE image a subset may correlate over -- issue #19, phase 2
	//of #14. A subset beside a hole is otherwise correlated over the hole too, and what
	//shows through a hole is background that does not move with the specimen, so the
	//answer is pulled toward no movement while the point still reports as solved.
	//Measured on synthetic images with an exact answer: points within 4 px of a 40 px
	//hole's edge were wrong by 0.69 px on average, correlating above 0.9.
	//
	//Built once per reference image from a Shape2D: a pixel the shape contains is valid.
	//A pixel outside the image is valid as far as this mask is concerned -- the solvers'
	//own bounds checks already decide what happens there, and this mask is about what
	//the specimen shows, not about where the image ends.
	//
	//⚑ The invalid region is then grown by `pull_in` pixels (a square, Chebyshev
	//neighbourhood), so the valid region is pulled back from every edge. Masking alone
	//is not enough: the reference gradient's 4th-order stencil reaches 2 px across an
	//edge and the target's bicubic B-spline further, so pixels right at the edge still
	//carry the hole's intensities. Measured beside the same hole: masked without pulling
	//in, 0.063 px; pulled in 2 px, 0.008 px, the far-field error.
	class ValidPixelMask2D
	{
	public:
		ValidPixelMask2D(const Shape2D& valid, int width, int height, int pull_in);

		int width, height;
		int pull_in;

		//after pull-in; true outside the image
		bool isValid(int x, int y) const;

		//valid pixels in the (2*radius_x + 1) x (2*radius_y + 1) subset centred on
		//(center_x, center_y), and that count as a share of the whole subset
		int validCount(int center_x, int center_y, int radius_x, int radius_y) const;
		float validShare(int center_x, int center_y, int radius_x, int radius_y) const;

		//the same subset as a 0/1 weight, rows by y; returns the valid count
		int fillWeight(int center_x, int center_y, int radius_x, int radius_y, Eigen::MatrixXf& weight) const;

		//invalid pixels in the whole image, after pull-in
		int invalidCount() const;

	private:
		Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic> valid; //height x width, 1 = valid
	};

}//namespace opencorr

#endif //_MASK_H_
