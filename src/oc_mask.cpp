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

#include <algorithm>
#include <string>

#include "oc_mask.h"

namespace opencorr
{
	namespace
	{
		//true where any of `invalid` lies within `reach` along one axis: a running count
		//over the window, so the cost does not grow with the reach. Pixels beyond the
		//ends count as valid.
		void growAlongRows(const Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& invalid,
			Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& grown, int reach)
		{
			const int rows = (int)invalid.rows(), cols = (int)invalid.cols();
			grown.setZero(rows, cols);
			for (int r = 0; r < rows; r++)
			{
				int count = 0;
				for (int c = 0; c < std::min(reach, cols); c++)
				{
					count += invalid(r, c);
				}
				for (int c = 0; c < cols; c++)
				{
					if (c + reach < cols) count += invalid(r, c + reach);
					if (c - reach - 1 >= 0) count -= invalid(r, c - reach - 1);
					grown(r, c) = count > 0 ? 1 : 0;
				}
			}
		}
	}

	ValidPixelMask2D::ValidPixelMask2D(const Shape2D& valid_region, int width, int height, int pull_in)
		: width(width), height(height), pull_in(pull_in)
	{
		if (width < 1 || height < 1 || pull_in < 0)
		{
			throw std::string("ValidPixelMask2D: the image must have pixels, and pull_in must not be negative");
		}

		Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic> invalid(height, width);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				invalid(y, x) = valid_region.contains(x, y) ? 0 : 1;
			}
		}

		//a square neighbourhood is separable: grow along rows, then along columns
		if (pull_in > 0)
		{
			Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic> along_rows, along_both;
			growAlongRows(invalid, along_rows, pull_in);
			Eigen::Array<uint8_t, Eigen::Dynamic, Eigen::Dynamic> transposed = along_rows.transpose();
			growAlongRows(transposed, along_both, pull_in);
			invalid = along_both.transpose();
		}

		valid = 1 - invalid;
	}

	bool ValidPixelMask2D::isValid(int x, int y) const
	{
		if (x < 0 || y < 0 || x >= width || y >= height)
		{
			return true;
		}
		return valid(y, x) != 0;
	}

	int ValidPixelMask2D::validCount(int center_x, int center_y, int radius_x, int radius_y) const
	{
		int count = 0;
		for (int y = center_y - radius_y; y <= center_y + radius_y; y++)
		{
			for (int x = center_x - radius_x; x <= center_x + radius_x; x++)
			{
				count += isValid(x, y) ? 1 : 0;
			}
		}
		return count;
	}

	float ValidPixelMask2D::validShare(int center_x, int center_y, int radius_x, int radius_y) const
	{
		const int size = (2 * radius_x + 1) * (2 * radius_y + 1);
		return (float)validCount(center_x, center_y, radius_x, radius_y) / (float)size;
	}

	int ValidPixelMask2D::fillWeight(int center_x, int center_y, int radius_x, int radius_y, Eigen::MatrixXf& weight) const
	{
		weight.resize(2 * radius_y + 1, 2 * radius_x + 1);
		int count = 0;
		for (int r = 0; r < weight.rows(); r++)
		{
			for (int c = 0; c < weight.cols(); c++)
			{
				const bool ok = isValid(center_x + c - radius_x, center_y + r - radius_y);
				weight(r, c) = ok ? 1.f : 0.f;
				count += ok ? 1 : 0;
			}
		}
		return count;
	}

	int ValidPixelMask2D::invalidCount() const
	{
		return (int)(valid == 0).count();
	}

}//namespace opencorr
