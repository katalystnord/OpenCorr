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
	Cine2D::Cine2D(const std::string& file_path)
	{
		//TO_8_BIT is only meaningful for 16-bit source cines (it sets the
		//linear scale factor down to 8 bits); 10-bit-packed sources need one
		//of the dedicated 10-bit conversion modes to land in the same 8-bit
		//range. Bit depth isn't known until the header is read, so probe it
		//with a throwaway header-only construction first, then reconstruct
		//with the right conversion mode if needed. 8-bit source cines are
		//unaffected by conversion_type either way.
		hc = std::make_unique<hypercine::HyperCine>(file_path.c_str(), hypercine::HyperCine::TO_8_BIT);
		if (hc->bit_depth() == hypercine::HyperCine::BIT_DEPTH_10_PACKED)
		{
			hc = std::make_unique<hypercine::HyperCine>(file_path.c_str(), hypercine::HyperCine::LINEAR_10_TO_8);
		}
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
		for (int r = 0; r < h; r++)
		{
			for (int c = 0; c < w; c++)
			{
				float value = frame_data[r * w + c];
				value = value < 0.f ? 0.f : (value > 255.f ? 255.f : value);
				image.cv_mat.at<uchar>(r, c) = (uchar)(value + 0.5f);
			}
		}
		cv::cv2eigen(image.cv_mat, image.eg_mat);
		image.file_path = hc->file_name() + "#" + std::to_string(frame_id);

		return image;
	}

}//namespace opencorr
