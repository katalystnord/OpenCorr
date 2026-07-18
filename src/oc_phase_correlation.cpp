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

#include "oc_phase_correlation.h"

namespace opencorr
{
	PhaseCorrelation2D::PhaseCorrelation2D(int width, int height)
		: width(width), height(height), hamming_window(true)
	{
		int size = width * height;
		buf_a = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * size);
		buf_b = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * size);
		buf_result = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * size);

		//in-place 2D complex-to-complex transforms; FFTW does the row+column pass
		//internally in one call rather than DICe's manual separable row-then-column loop
		plan_forward_a = fftwf_plan_dft_2d(height, width, buf_a, buf_a, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_forward_b = fftwf_plan_dft_2d(height, width, buf_b, buf_b, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_inverse = fftwf_plan_dft_2d(height, width, buf_result, buf_result, FFTW_BACKWARD, FFTW_ESTIMATE);
	}

	PhaseCorrelation2D::~PhaseCorrelation2D()
	{
		fftwf_destroy_plan(plan_forward_a);
		fftwf_destroy_plan(plan_forward_b);
		fftwf_destroy_plan(plan_inverse);
		fftwf_free(buf_a);
		fftwf_free(buf_b);
		fftwf_free(buf_result);
	}

	bool PhaseCorrelation2D::getHammingWindow() const
	{
		return hamming_window;
	}

	void PhaseCorrelation2D::setHammingWindow(bool apply_hamming_window)
	{
		hamming_window = apply_hamming_window;
	}

	void PhaseCorrelation2D::fillWindowed(Image2D& image, fftwf_complex* buf)
	{
		//avoid M_PI: MSVC's <cmath> only defines it when _USE_MATH_DEFINES is set
		//before the include, which nothing in this build sets (same fix as oc_uncertainty.cpp)
		static const float two_pi = 2.f * (float)std::acos(-1.0);

		for (int y = 0; y < height; y++)
		{
			float wy = hamming_window ? (0.54f - 0.46f * std::cos(two_pi * y / (height - 1))) : 1.f;
			for (int x = 0; x < width; x++)
			{
				float wx = hamming_window ? (0.54f - 0.46f * std::cos(two_pi * x / (width - 1))) : 1.f;
				buf[y * width + x][0] = image.eg_mat(y, x) * wx * wy;
				buf[y * width + x][1] = 0.f;
			}
		}
	}

	float PhaseCorrelation2D::compute(Image2D& ref_img, Image2D& tar_img, float& u, float& v)
	{
		fillWindowed(ref_img, buf_a);
		fillWindowed(tar_img, buf_b);

		fftwf_execute(plan_forward_a);
		fftwf_execute(plan_forward_b);

		//cross-power spectrum, normalized by its own magnitude -- this normalization is
		//what turns plain cross-correlation into phase correlation: it whitens the
		//spectrum so the result depends only on phase (relative shift) rather than
		//magnitude (local contrast/texture strength), giving a sharper, more reliable
		//peak for large-motion/low-texture cases than plain cross-correlation would
		int size = width * height;
		for (int i = 0; i < size; i++)
		{
			//cross = A * conj(B)
			float cross_r = buf_a[i][0] * buf_b[i][0] + buf_a[i][1] * buf_b[i][1];
			float cross_i = buf_a[i][1] * buf_b[i][0] - buf_a[i][0] * buf_b[i][1];
			float mag = std::sqrt(cross_r * cross_r + cross_i * cross_i);
			if (mag > 1e-10f)
			{
				buf_result[i][0] = cross_r / mag;
				buf_result[i][1] = cross_i / mag;
			}
			else
			{
				buf_result[i][0] = 0.f;
				buf_result[i][1] = 0.f;
			}
		}

		fftwf_execute(plan_inverse);
		//FFTW's inverse transform is unnormalized (a forward+inverse round trip scales
		//values by width*height); normalize so the returned peak magnitude is a stable,
		//comparable confidence score independent of image size
		float norm = 1.f / (float)size;

		//find the peak magnitude location; skip a spurious peak sitting exactly at (0,0)
		//if a comparable peak exists elsewhere -- ported from DICe's own workaround for a
		//false zero-shift peak this normalization can otherwise produce
		float max_real = 0.f, next_real = 0.f;
		int peak_x = 0, peak_y = 0, next_x = 0, next_y = 0;
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				float val = std::fabs(buf_result[y * width + x][0]) * norm;
				if (val > max_real)
				{
					max_real = val;
					peak_x = x;
					peak_y = y;
				}
				if (x == 0 && y == 0) continue;
				if (val > next_real)
				{
					next_real = val;
					next_x = x;
					next_y = y;
				}
			}
		}
		if (peak_x != next_x && next_x > 1)
		{
			peak_x = next_x;
			max_real = next_real;
		}
		if (peak_y != next_y && next_y > 1)
		{
			peak_y = next_y;
			max_real = next_real;
		}

		//convert wrapped FFT bin index to a signed pixel displacement
		u = (peak_x >= width / 2) ? (float)(width - peak_x) : -(float)peak_x;
		v = (peak_y >= height / 2) ? (float)(height - peak_y) : -(float)peak_y;

		return max_real;
	}

}//namespace opencorr
