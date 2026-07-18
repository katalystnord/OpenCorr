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

#ifndef _CINE_H_
#define _CINE_H_

#include <memory>
#include <string>

#include "oc_image.h"

namespace hypercine { class HyperCine; }

namespace opencorr
{
	//Reads .cine files (Phantom/Vision Research high-speed camera format) as
	//a sequence of Image2D frames, via a vendored copy of hypercine
	//(dicengine/hypercine, BSD-3-Clause, deps/hypercine/). OpenCorr has no
	//native image-sequence abstraction (each example constructs one Image2D
	//per file path), so this class is scoped narrowly to file decoding only
	//-- looping over a sequence of frames remains the caller's job, same as
	//the existing convention for ordinary image files.
	//
	//16-bit source cines are normalized to an 8-bit intensity range during
	//decode (hypercine's TO_8_BIT conversion mode), matching Image2D's
	//existing 8-bit grayscale convention (Image2D loads every format via
	//cv::IMREAD_GRAYSCALE); 8-bit sources pass through unchanged. hypercine
	//also handles the 16-bit format's upside-down row storage internally,
	//so frames come out right-side-up.
	//
	//10-bit-packed sources need one of two different nonlinear-to-8-bit
	//lookup curves (see Cine10BitConversion below); which one is correct
	//depends on the recording, and there's no way to detect that from the
	//file itself, so it's a caller-supplied choice rather than a guess baked
	//in silently. The default (Linear) is the one verified against this
	//repo's own bundled 10-bit-packed test fixture (examples/cine/
	//packed_12bpp.cine, see cine_smoke_test.cpp) -- if a real-world 10-bit
	//recording's decoded intensities look visibly wrong (e.g. an
	//unexpectedly compressed or expanded contrast range), try Quadratic.

	enum class Cine10BitConversion
	{
		Linear, //hypercine's LINEAR_10_TO_8
		Quadratic //hypercine's QUAD_10_TO_8 (the same lookup table TO_8_BIT also uses for 10-bit-packed sources)
	};

	class Cine2D
	{
	private:
		std::unique_ptr<hypercine::HyperCine> hc;

	public:
		explicit Cine2D(const std::string& file_path, Cine10BitConversion ten_bit_conversion = Cine10BitConversion::Linear);
		~Cine2D();

		int width() const;
		int height() const;
		int firstFrameId() const; //cine frame ids are not necessarily 0-based
		int frameCount() const;

		Image2D getFrame(int frame_id); //decode a single frame as an 8-bit Image2D
	};

}//namespace opencorr

#endif //_CINE_H_
