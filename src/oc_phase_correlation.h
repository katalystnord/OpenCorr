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

#ifndef _PHASE_CORRELATION_H_
#define _PHASE_CORRELATION_H_

#include "fftw3.h"

#include "oc_image.h"

namespace opencorr
{
	//Whole-image FFT phase correlation: a coarse, large-motion translation
	//estimate between two full images, ported from DICe::phase_correlate_x_y()
	//(dicengine/dice, src/fft/DICe_FFT.cpp, BSD-3-Clause). Algorithm is a
	//direct port (Hamming-windowed FFT of both images, cross-power spectrum
	//normalized by its own magnitude -- the normalization step that makes
	//this phase correlation rather than plain cross-correlation -- inverse
	//FFT, peak location with signed-wraparound conversion, plus DICe's
	//spurious-(0,0)-peak workaround).
	//
	//Deliberately uses OpenCorr's existing FFTW dependency (already linked
	//in for FFTCC2D, see oc_fftcc.h) instead of vendoring DICe's bundled
	//KissFFT: FFTW does the 2D transform natively in one call rather than
	//DICe's manual separable row-then-column loop over 1D KissFFT calls, and
	//it avoids carrying two different FFT libraries for two conceptually
	//similar features.
	//
	//Complements, not replaces, OpenCorr's existing initial-guess mechanisms
	//(FFTCC2D, FeatureAffine2D): those work at the scale of one subset window;
	//this works at the scale of the whole image, which is useful specifically
	//when there isn't enough local texture in an individual subset for
	//FFTCC2D/FeatureAffine2D to find a reliable local match (e.g. large
	//rigid-body motion on a low-texture sample), but the image as a whole
	//still has enough content to phase-correlate.
	//
	//Not thread-safe: unlike OpenCorr's per-POI solvers (ICGN2D1 etc., which pool
	//per-thread scratch state because they're called from inside an OpenMP loop over
	//POIs), PhaseCorrelation2D is meant to run once per image pair as a whole-image
	//initializer BEFORE the per-POI parallel stage, not inside it -- so it keeps a
	//single set of FFTW buffers/plans as instance state rather than pooling them. If a
	//caller does need concurrent calls (e.g. processing several independent image pairs
	//in parallel), construct one instance per thread.

	class PhaseCorrelation2D
	{
	public:
		PhaseCorrelation2D(int width, int height);
		~PhaseCorrelation2D();

		bool getHammingWindow() const;
		void setHammingWindow(bool apply_hamming_window); //default true, reduces edge-induced spectral leakage

		//estimates the whole-image translation (u, v) that best phase-aligns tar_img to
		//ref_img (both must be width x height, matching the constructor); returns the
		//cross-power-spectrum peak magnitude as a rough confidence score -- low values
		//(near zero) mean no clear dominant translation was found (e.g. the images are
		//identical, or too dissimilar/low-texture for phase correlation to lock onto
		//anything)
		float compute(Image2D& ref_img, Image2D& tar_img, float& u, float& v);

	private:
		int width, height;
		bool hamming_window;

		fftwf_complex* buf_a;
		fftwf_complex* buf_b;
		fftwf_complex* buf_result;
		fftwf_plan plan_forward_a;
		fftwf_plan plan_forward_b;
		fftwf_plan plan_inverse;

		void fillWindowed(Image2D& image, fftwf_complex* buf);
	};

}//namespace opencorr

#endif //_PHASE_CORRELATION_H_
