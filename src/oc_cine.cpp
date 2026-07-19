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

#ifndef USE_FLOAT_STORAGE
#define USE_FLOAT_STORAGE 1 //hypercine's storage_t; matches Image2D's Eigen::MatrixXf (float) representation
#endif
#include "hypercine.h"

#include "oc_cine.h"

namespace opencorr
{
	Cine2D::Cine2D(const std::string& file_path, Cine10BitConversion ten_bit_conversion)
	{
		//TO_8_BIT and 10-bit-packed sources: contrary to what an earlier version of this
		//comment claimed, TO_8_BIT is NOT 16-bit-only -- hypercine.cpp maps it to the same
		//quad_10bit_to_8bit lookup table QUAD_10_TO_8 uses for 10-bit-packed sources too
		//(hypercine.cpp:353-354). So for 10-bit-packed files, TO_8_BIT and Quadratic are the
		//same curve; Linear is a genuinely different one (LINEAR_10_TO_8, its own lookup
		//table). Which curve is actually correct depends on the recording and isn't
		//knowable from the file itself -- see Cine10BitConversion's doc comment in oc_cine.h.
		//
		//Bit depth isn't known until the header is read, so probe it with a throwaway
		//header-only construction first (using TO_8_BIT, harmless for 8-bit/16-bit sources
		//and equivalent to Quadratic if the file turns out to be 10-bit-packed), then
		//reconstruct with the caller's requested curve if it turns out to be 10-bit-packed.
		hc = std::make_unique<hypercine::HyperCine>(file_path.c_str(), hypercine::HyperCine::TO_8_BIT);
		if (hc->bit_depth() == hypercine::HyperCine::BIT_DEPTH_10_PACKED
			&& ten_bit_conversion == Cine10BitConversion::Linear)
		{
			hc = std::make_unique<hypercine::HyperCine>(file_path.c_str(), hypercine::HyperCine::LINEAR_10_TO_8);
		}
		//else: already open with the quadratic curve (TO_8_BIT), which is what Quadratic
		//requests too -- no second construction needed
	}

	Cine2D::~Cine2D() = default;

	int Cine2D::width() const
	{
		return hc->width();
	}

	int Cine2D::height() const
	{
		return hc->height();
	}

	int Cine2D::firstFrameId() const
	{
		return hc->file_first_frame_id();
	}

	int Cine2D::frameCount() const
	{
		return hc->file_frame_count();
	}

	Image2D Cine2D::getFrame(int frame_id)
	{
		int w = hc->width();
		int h = hc->height();
		std::vector<storage_t> frame_data = hc->get_frame(frame_id); //storage_t = float, already normalized to [0, 255]-ish range

		Image2D image(w, h);
		//non-owning view over frame_data's own buffer (valid for the lifetime of this
		//call, which is all convertTo() needs) -- one vectorized OpenCV call instead of a
		//per-pixel loop. convertTo()'s saturate_cast<uchar> rounds to nearest and clamps
		//to [0, 255] internally, exactly matching the manual clamp-then-round this
		//replaces. OPENCV_DATA_TYPE (hypercine.h) is defined alongside storage_t's own
		//typedef to always name the matching OpenCV type (double/float/uint16_t storage
		//select USE_DOUBLE_STORAGE/USE_FLOAT_STORAGE/USE_INT_STORAGE at build time --
		//OpenCorr's own CMakeLists.txt sets USE_FLOAT_STORAGE), so this stays correct
		//regardless of which one is actually configured.
		cv::Mat frame_mat(h, w, OPENCV_DATA_TYPE, frame_data.data());
		frame_mat.convertTo(image.cv_mat, CV_8U);
		cv::cv2eigen(image.cv_mat, image.eg_mat);
		image.file_path = hc->file_name() + "#" + std::to_string(frame_id);

		return image;
	}

}//namespace opencorr
