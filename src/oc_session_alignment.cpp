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
#include <cmath>
#include <memory>

#include <opencv2/calib3d.hpp>

#include "oc_icgn.h"
#include "oc_mask.h"
#include "oc_session_alignment.h"
#include "oc_sift.h"

namespace opencorr
{
	Point2D Similarity2D::map(const Point2D& reference) const
	{
		const double c = std::cos(angle), s = std::sin(angle);
		return Point2D((float)(scale * (c * reference.x - s * reference.y) + tx),
			(float)(scale * (s * reference.x + c * reference.y) + ty));
	}

	Point2D Similarity2D::unmap(const Point2D& target) const
	{
		const double c = std::cos(angle), s = std::sin(angle);
		const double rx = (target.x - tx) / scale, ry = (target.y - ty) / scale;
		return Point2D((float)(c * rx + s * ry), (float)(-s * rx + c * ry));
	}

	void Similarity2D::seed(POI2D& poi) const
	{
		const double c = std::cos(angle), s = std::sin(angle);
		const Point2D x = map(Point2D(poi.x, poi.y));
		for (float& p : poi.deformation.p) p = 0.f;
		poi.deformation.u = x.x - poi.x;
		poi.deformation.v = x.y - poi.y;
		poi.deformation.ux = (float)(scale * c - 1.0);
		poi.deformation.uy = (float)(-scale * s);
		poi.deformation.vx = (float)(scale * s);
		poi.deformation.vy = (float)(scale * c - 1.0);
	}

	void Similarity2D::toReference(POI2D& poi) const
	{
		const double c = std::cos(angle), s = std::sin(angle);

		//the material point's position in the target, back through the camera
		const Point2D at = unmap(Point2D(poi.x + poi.deformation.u, poi.y + poi.deformation.v));
		poi.deformation.u = at.x - poi.x;
		poi.deformation.v = at.y - poi.y;

		//G = I + grad u maps reference to target; the specimen's own is R^T G / scale
		const double g11 = 1.0 + poi.deformation.ux, g12 = poi.deformation.uy;
		const double g21 = poi.deformation.vx, g22 = 1.0 + poi.deformation.vy;
		poi.deformation.ux = (float)((c * g11 + s * g21) / scale - 1.0);
		poi.deformation.uy = (float)((c * g12 + s * g22) / scale);
		poi.deformation.vx = (float)((-s * g11 + c * g21) / scale);
		poi.deformation.vy = (float)((-s * g12 + c * g22) / scale - 1.0);
	}

	SessionAlignment2D::SessionAlignment2D(int subset_radius, int grid_step, int thread_number)
		: subset_radius(subset_radius), grid_step(grid_step), thread_number(thread_number) {}

	void SessionAlignment2D::setImages(Image2D& ref_img, Image2D& tar_img)
	{
		this->ref_img = &ref_img;
		this->tar_img = &tar_img;
	}

	void SessionAlignment2D::setGround(const Shape2D* ground) { this->ground = ground; }

	void SessionAlignment2D::setMinimums(int min_matches, int min_refined)
	{
		this->min_matches = min_matches;
		this->min_refined = min_refined;
	}

	namespace
	{
		//least-squares similarity through point pairs: x = a X - b Y + tx, y = b X + a Y + ty,
		//about the centroid so the normal equations stay well conditioned
		bool fitSimilarity(const std::vector<Point2D>& from, const std::vector<Point2D>& to,
			const std::vector<char>& keep, Similarity2D& out)
		{
			double fx = 0, fy = 0, tx = 0, ty = 0;
			int n = 0;
			for (size_t i = 0; i < from.size(); i++)
			{
				if (!keep[i]) continue;
				fx += from[i].x; fy += from[i].y; tx += to[i].x; ty += to[i].y; n++;
			}
			if (n < 2) return false;
			fx /= n; fy /= n; tx /= n; ty /= n;

			double sxx = 0, dot = 0, cross = 0;
			for (size_t i = 0; i < from.size(); i++)
			{
				if (!keep[i]) continue;
				const double X = from[i].x - fx, Y = from[i].y - fy;
				const double x = to[i].x - tx, y = to[i].y - ty;
				sxx += X * X + Y * Y;
				dot += X * x + Y * y;
				cross += X * y - Y * x;
			}
			if (sxx <= 0) return false;
			const double a = dot / sxx, b = cross / sxx;
			out.scale = std::hypot(a, b);
			out.angle = std::atan2(b, a);
			out.tx = tx - (a * fx - b * fy);
			out.ty = ty - (b * fx + a * fy);
			return out.scale > 0;
		}

		double residual(const Similarity2D& s, const Point2D& from, const Point2D& to)
		{
			const Point2D m = s.map(from);
			return std::hypot(m.x - to.x, m.y - to.y);
		}
	}

	AlignmentResult SessionAlignment2D::compute()
	{
		AlignmentResult result;
		if (ref_img == nullptr || tar_img == nullptr) return result;
		result.status = AlignmentStatus::too_few_matches;
		if (ref_img->cv_mat.empty() || tar_img->cv_mat.empty()) return result;

		//1. SIFT matches on the ground
		SIFT2D sift;
		sift.setImages(*ref_img, *tar_img);
		sift.prepare();
		sift.compute();

		std::vector<cv::Point2f> from, to;
		for (size_t i = 0; i < sift.ref_matched_kp.size(); i++)
		{
			const Point2D& r = sift.ref_matched_kp[i];
			if (ground && !ground->contains((int)std::lround(r.x), (int)std::lround(r.y))) continue;
			from.emplace_back(r.x, r.y);
			to.emplace_back(sift.tar_matched_kp[i].x, sift.tar_matched_kp[i].y);
		}
		result.matches = (int)from.size();
		if (result.matches < std::max(min_matches, 2)) return result;

		//2. RANSAC similarity
		cv::Mat inliers;
		const cv::Mat h = cv::estimateAffinePartial2D(from, to, inliers, cv::RANSAC, 1.0);
		if (h.empty()) return result;
		result.inliers = cv::countNonZero(inliers);
		if (result.inliers < std::max(min_matches, 2)) return result;

		Similarity2D coarse;
		coarse.scale = std::hypot(h.at<double>(0, 0), h.at<double>(1, 0));
		coarse.angle = std::atan2(h.at<double>(1, 0), h.at<double>(0, 0));
		coarse.tx = h.at<double>(0, 2);
		coarse.ty = h.at<double>(1, 2);

		//3. refinement on a grid of subsets lying wholly on the ground
		const int r = subset_radius, side = 2 * r + 1;
		std::unique_ptr<ValidPixelMask2D> on_ground;
		if (ground) on_ground = std::make_unique<ValidPixelMask2D>(*ground, ref_img->width, ref_img->height, 0);

		std::vector<POI2D> queue;
		for (int y = r + 1; y < ref_img->height - r - 1; y += grid_step)
			for (int x = r + 1; x < ref_img->width - r - 1; x += grid_step)
			{
				if (on_ground && on_ground->validCount(x, y, r, r) < side * side) continue;
				POI2D poi(x, y);
				coarse.seed(poi);
				const Point2D t(x + poi.deformation.u, y + poi.deformation.v);
				if (t.x < r + 1 || t.y < r + 1 || t.x > tar_img->width - r - 2 || t.y > tar_img->height - r - 2) continue;
				queue.push_back(poi);
			}
		result.refined_points = (int)queue.size();

		result.status = AlignmentStatus::too_few_ground_points;
		if (result.refined_points < min_refined) return result;

		ICGN2D1 icgn(r, r, 0.001f, 20, thread_number);
		icgn.setImages(*ref_img, *tar_img);
		icgn.prepare();
		icgn.compute(queue);

		std::vector<Point2D> centre, found;
		for (const auto& poi : queue)
		{
			if (!(poi.result.zncc >= 0.9f) || !std::isfinite(poi.deformation.u) || !std::isfinite(poi.deformation.v)) continue;
			centre.push_back(Point2D(poi.x, poi.y));
			found.push_back(Point2D(poi.x + poi.deformation.u, poi.y + poi.deformation.v));
		}
		if ((int)centre.size() < min_refined) return result;

		//reweighted: keep what the current fit describes, refit, until the kept set
		//settles. The threshold is three times the kept points' median residual, with a
		//floor of 0.02 px so a perfect fit does not start rejecting its own noise.
		Similarity2D fine = coarse;
		std::vector<char> keep(centre.size(), 1);
		for (int round = 0; round < 20; round++)
		{
			std::vector<double> res(centre.size());
			for (size_t i = 0; i < centre.size(); i++) res[i] = residual(fine, centre[i], found[i]);
			std::vector<double> kept_res;
			for (size_t i = 0; i < centre.size(); i++) if (keep[i]) kept_res.push_back(res[i]);
			if (kept_res.empty()) return result;
			std::nth_element(kept_res.begin(), kept_res.begin() + kept_res.size() / 2, kept_res.end());
			const double limit = std::max(3.0 * kept_res[kept_res.size() / 2], 0.02);

			std::vector<char> next(centre.size());
			for (size_t i = 0; i < centre.size(); i++) next[i] = res[i] <= limit;
			const bool settled = round > 0 && next == keep;
			keep = next;
			if (settled) break;
			if (!fitSimilarity(centre, found, keep, fine)) return result;
		}

		double ss = 0;
		int used = 0;
		for (size_t i = 0; i < centre.size(); i++)
		{
			if (!keep[i]) continue;
			const double e = residual(fine, centre[i], found[i]);
			ss += e * e;
			used++;
		}
		result.refined_used = used;
		if (used < min_refined) return result;
		result.residual_rms = (float)std::sqrt(ss / used);
		result.transform = fine;
		result.status = AlignmentStatus::ok;
		return result;
	}

}//namespace opencorr
