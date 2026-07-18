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
	//10-bit-packed and 16-bit source cines are both normalized to an 8-bit
	//intensity range during decode (hypercine's TO_8_BIT conversion mode),
	//matching Image2D's existing 8-bit grayscale convention (Image2D loads
	//every format via cv::IMREAD_GRAYSCALE); 8-bit source cines pass through
	//unchanged. hypercine also handles the 16-bit format's upside-down row
	//storage internally, so frames come out right-side-up.

	class Cine2D
	{
	private:
		std::unique_ptr<hypercine::HyperCine> hc;

	public:
		explicit Cine2D(const std::string& file_path);
		~Cine2D();

		int width() const;
		int height() const;
		int firstFrameId() const; //cine frame ids are not necessarily 0-based
		int frameCount() const;

		Image2D getFrame(int frame_id); //decode a single frame as an 8-bit Image2D
	};

}//namespace opencorr

#endif //_CINE_H_
