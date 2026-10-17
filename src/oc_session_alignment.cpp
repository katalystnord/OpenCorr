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
#include <limits>
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

	Homography2D Homography2D::from(const Similarity2D& similarity)
	{
		const double c = similarity.scale * std::cos(similarity.angle), s = similarity.scale * std::sin(similarity.angle);
		Homography2D out;
		const double h[9] = { c, -s, similarity.tx, s, c, similarity.ty, 0, 0, 1 };
		for (int i = 0; i < 9; i++) out.h[i] = h[i];
		return out;
	}

	Point2D Homography2D::map(const Point2D& reference) const
	{
		const double X = reference.x, Y = reference.y;
		const double w = h[6] * X + h[7] * Y + h[8];
		return Point2D((float)((h[0] * X + h[1] * Y + h[2]) / w), (float)((h[3] * X + h[4] * Y + h[5]) / w));
	}

	Point2D Homography2D::unmap(const Point2D& target) const
	{
		//the adjugate: the inverse up to a scale the division removes
		const double a = h[4] * h[8] - h[5] * h[7], b = h[2] * h[7] - h[1] * h[8], c = h[1] * h[5] - h[2] * h[4];
		const double d = h[5] * h[6] - h[3] * h[8], e = h[0] * h[8] - h[2] * h[6], f = h[2] * h[3] - h[0] * h[5];
		const double g = h[3] * h[7] - h[4] * h[6], k = h[1] * h[6] - h[0] * h[7], l = h[0] * h[4] - h[1] * h[3];
		const double x = target.x, y = target.y;
		const double w = g * x + k * y + l;
		return Point2D((float)((a * x + b * y + c) / w), (float)((d * x + e * y + f) / w));
	}

	void Homography2D::jacobian(const Point2D& reference, double j[4]) const
	{
		const double X = reference.x, Y = reference.y;
		const double w = h[6] * X + h[7] * Y + h[8];
		const double x = (h[0] * X + h[1] * Y + h[2]) / w, y = (h[3] * X + h[4] * Y + h[5]) / w;
		j[0] = (h[0] - x * h[6]) / w;
		j[1] = (h[1] - x * h[7]) / w;
		j[2] = (h[3] - y * h[6]) / w;
		j[3] = (h[4] - y * h[7]) / w;
	}

	void Homography2D::seed(POI2D& poi) const
	{
		const Point2D x = map(Point2D(poi.x, poi.y));
		double j[4];
		jacobian(Point2D(poi.x, poi.y), j);
		for (float& p : poi.deformation.p) p = 0.f;
		poi.deformation.u = x.x - poi.x;
		poi.deformation.v = x.y - poi.y;
		poi.deformation.ux = (float)(j[0] - 1.0);
		poi.deformation.uy = (float)j[1];
		poi.deformation.vx = (float)j[2];
		poi.deformation.vy = (float)(j[3] - 1.0);
	}

	void Homography2D::toReference(POI2D& poi) const
	{
		//the material point's position in the target, back through the camera
		const Point2D at = unmap(Point2D(poi.x + poi.deformation.u, poi.y + poi.deformation.v));

		//G = I + grad u maps reference to target, and is the camera's local J at where the
		//material point now sits, times the specimen's own F: F = J(at)^-1 G
		double j[4];
		jacobian(at, j);
		const double det = j[0] * j[3] - j[1] * j[2];
		const double i11 = j[3] / det, i12 = -j[1] / det, i21 = -j[2] / det, i22 = j[0] / det;
		const double g11 = 1.0 + poi.deformation.ux, g12 = poi.deformation.uy;
		const double g21 = poi.deformation.vx, g22 = 1.0 + poi.deformation.vy;

		poi.deformation.u = at.x - poi.x;
		poi.deformation.v = at.y - poi.y;
		poi.deformation.ux = (float)(i11 * g11 + i12 * g21 - 1.0);
		poi.deformation.uy = (float)(i11 * g12 + i12 * g22);
		poi.deformation.vx = (float)(i21 * g11 + i22 * g21);
		poi.deformation.vy = (float)(i21 * g12 + i22 * g22 - 1.0);
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

	void SessionAlignment2D::setKnownLength(const Point2D& a, const Point2D& b)
	{
		known_length[0] = Point2D((float)std::lround(a.x), (float)std::lround(a.y));
		known_length[1] = Point2D((float)std::lround(b.x), (float)std::lround(b.y));
		has_known_length = true;
	}

	void SessionAlignment2D::clearKnownLength() { has_known_length = false; }

	void SessionAlignment2D::setPerspective(bool on) { perspective = on; }

	namespace
	{
		//least-squares similarity through point pairs: x = a X - b Y + tx, y = b X + a Y + ty,
		//about the centroid so the normal equations stay well conditioned
		//fixed_scale above zero holds the scale at that value and fits only the rotation
		//and the shift, which is the least-squares rotation of the centred pairs
		bool fitSimilarity(const std::vector<Point2D>& from, const std::vector<Point2D>& to,
			const std::vector<char>& keep, Similarity2D& out, double fixed_scale = 0.0)
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
			double a = dot / sxx, b = cross / sxx;
			if (fixed_scale > 0)
			{
				const double turn = std::atan2(cross, dot);
				a = fixed_scale * std::cos(turn);
				b = fixed_scale * std::sin(turn);
			}
			out.scale = std::hypot(a, b);
			out.angle = std::atan2(b, a);
			out.tx = tx - (a * fx - b * fy);
			out.ty = ty - (b * fx + a * fy);
			return out.scale > 0;
		}

		template <class Camera>
		double residual(const Camera& s, const Point2D& from, const Point2D& to)
		{
			const Point2D m = s.map(from);
			return std::hypot(m.x - to.x, m.y - to.y);
		}

		//least-squares homography through point pairs, minimising the distance in the
		//target (the geometric error, not the algebraic one), by Gauss-Newton from the
		//homography given. About the centroids, scaled to unit spread, so the eight
		//parameters are of one size and the normal equations stay well conditioned.
		bool fitHomography(const std::vector<Point2D>& from, const std::vector<Point2D>& to,
			const std::vector<char>& keep, Homography2D& out)
		{
			double fx = 0, fy = 0, tx = 0, ty = 0;
			int n = 0;
			for (size_t i = 0; i < from.size(); i++)
			{
				if (!keep[i]) continue;
				fx += from[i].x; fy += from[i].y; tx += to[i].x; ty += to[i].y; n++;
			}
			if (n < 4) return false;
			fx /= n; fy /= n; tx /= n; ty /= n;
			double fs = 0, ts = 0;
			for (size_t i = 0; i < from.size(); i++)
			{
				if (!keep[i]) continue;
				fs += std::hypot(from[i].x - fx, from[i].y - fy);
				ts += std::hypot(to[i].x - tx, to[i].y - ty);
			}
			fs = fs / n; ts = ts / n;
			if (!(fs > 0) || !(ts > 0)) return false;

			//normalised: Hn = Tt H Tf^-1, with Tf, Tt centring and scaling
			const Eigen::Matrix3d Tf = (Eigen::Matrix3d() << 1 / fs, 0, -fx / fs, 0, 1 / fs, -fy / fs, 0, 0, 1).finished();
			const Eigen::Matrix3d Tt = (Eigen::Matrix3d() << 1 / ts, 0, -tx / ts, 0, 1 / ts, -ty / ts, 0, 0, 1).finished();
			Eigen::Matrix3d H;
			for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) H(r, c) = out.h[3 * r + c];
			Eigen::Matrix3d Hn = Tt * H * Tf.inverse();
			Hn /= Hn(2, 2);

			std::vector<double> X, Y, x, y;
			for (size_t i = 0; i < from.size(); i++)
			{
				if (!keep[i]) continue;
				X.push_back((from[i].x - fx) / fs); Y.push_back((from[i].y - fy) / fs);
				x.push_back((to[i].x - tx) / ts); y.push_back((to[i].y - ty) / ts);
			}

			for (int iteration = 0; iteration < 10; iteration++)
			{
				Eigen::MatrixXd J(2 * n, 8);
				Eigen::VectorXd r(2 * n);
				for (int i = 0; i < n; i++)
				{
					const double w = Hn(2, 0) * X[i] + Hn(2, 1) * Y[i] + 1;
					const double px = (Hn(0, 0) * X[i] + Hn(0, 1) * Y[i] + Hn(0, 2)) / w;
					const double py = (Hn(1, 0) * X[i] + Hn(1, 1) * Y[i] + Hn(1, 2)) / w;
					r(2 * i) = x[i] - px;
					r(2 * i + 1) = y[i] - py;
					J.row(2 * i) << X[i] / w, Y[i] / w, 1 / w, 0, 0, 0, -px * X[i] / w, -px * Y[i] / w;
					J.row(2 * i + 1) << 0, 0, 0, X[i] / w, Y[i] / w, 1 / w, -py * X[i] / w, -py * Y[i] / w;
				}
				const Eigen::VectorXd step = J.colPivHouseholderQr().solve(r);
				if (!step.allFinite()) return false;
				Hn(0, 0) += step(0); Hn(0, 1) += step(1); Hn(0, 2) += step(2);
				Hn(1, 0) += step(3); Hn(1, 1) += step(4); Hn(1, 2) += step(5);
				Hn(2, 0) += step(6); Hn(2, 1) += step(7);
				if (step.norm() < 1e-12) break;
			}

			H = Tt.inverse() * Hn * Tf;
			if (!H.allFinite() || H(2, 2) == 0) return false;
			H /= H(2, 2);
			for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) out.h[3 * r + c] = H(r, c);
			return true;
		}

		//reweighted: keep what the current fit describes, refit, until the kept set
		//settles. The threshold is three times the kept points' median residual, with a
		//floor of 0.02 px so a perfect fit does not start rejecting its own noise.
		template <class Camera, class Fit>
		bool reweighted(const std::vector<Point2D>& centre, const std::vector<Point2D>& found,
			std::vector<char>& keep, Camera& fine, Fit fit)
		{
			for (int round = 0; round < 20; round++)
			{
				std::vector<double> res(centre.size());
				for (size_t i = 0; i < centre.size(); i++) res[i] = residual(fine, centre[i], found[i]);
				std::vector<double> kept_res;
				for (size_t i = 0; i < centre.size(); i++) if (keep[i]) kept_res.push_back(res[i]);
				if (kept_res.empty()) return false;
				std::nth_element(kept_res.begin(), kept_res.begin() + kept_res.size() / 2, kept_res.end());
				const double limit = std::max(3.0 * kept_res[kept_res.size() / 2], 0.02);

				std::vector<char> next(centre.size());
				for (size_t i = 0; i < centre.size(); i++) next[i] = res[i] <= limit;
				const bool settled = round > 0 && next == keep;
				keep = next;
				if (settled) break;
				if (!fit(keep, fine)) return false;
			}
			return true;
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

		//asked for perspective, the refinement starts from a RANSAC homography: under a
		//tilt the similarity is wrong by pixels at the edges of the picture, beyond where
		//the subsets there would converge from
		Homography2D coarse_h = Homography2D::from(coarse);
		if (perspective)
		{
			const cv::Mat hh = cv::findHomography(from, to, cv::RANSAC, 1.0);
			if (!hh.empty())
				for (int i = 0; i < 9; i++) coarse_h.h[i] = hh.at<double>(i / 3, i % 3) / hh.at<double>(2, 2);
		}

		const int r = subset_radius, side = 2 * r + 1;
		ICGN2D1 icgn(r, r, 0.001f, 20, thread_number);
		icgn.setImages(*ref_img, *tar_img);
		icgn.prepare();

		auto fitsBoth = [&](const POI2D& poi) {
			const Point2D t(poi.x + poi.deformation.u, poi.y + poi.deformation.v);
			return poi.x >= r + 1 && poi.y >= r + 1 && poi.x <= ref_img->width - r - 2
				&& poi.y <= ref_img->height - r - 2 && t.x >= r + 1 && t.y >= r + 1
				&& t.x <= tar_img->width - r - 2 && t.y <= tar_img->height - r - 2;
		};

		//the known length: both ends correlated on the target, and the scale taken from them
		double fixed_scale = 0.0;
		if (has_known_length && !perspective)
		{
			result.status = AlignmentStatus::known_length_unmeasured;
			std::vector<POI2D> ends;
			for (const Point2D& end : known_length)
			{
				POI2D poi(end);
				coarse.seed(poi);
				if (!fitsBoth(poi)) return result;
				ends.push_back(poi);
			}
			icgn.compute(ends);
			for (int i = 0; i < 2; i++)
			{
				const POI2D& poi = ends[i];
				result.known_length_zncc[i] = poi.result.zncc;
				if (!(poi.result.zncc >= 0.9f) || !std::isfinite(poi.deformation.u) || !std::isfinite(poi.deformation.v))
					return result;
				result.known_length_target[i] = Point2D(poi.x + poi.deformation.u, poi.y + poi.deformation.v);
			}
			const double before = std::hypot(known_length[1].x - known_length[0].x, known_length[1].y - known_length[0].y);
			const double after = std::hypot(result.known_length_target[1].x - result.known_length_target[0].x,
				result.known_length_target[1].y - result.known_length_target[0].y);
			if (!(before > 0) || !(after > 0)) return result;
			fixed_scale = after / before;
			result.known_length_used = true;
		}

		//3. refinement on a grid of subsets lying wholly on the ground
		std::unique_ptr<ValidPixelMask2D> on_ground;
		if (ground) on_ground = std::make_unique<ValidPixelMask2D>(*ground, ref_img->width, ref_img->height, 0);

		std::vector<POI2D> queue;
		for (int y = r + 1; y < ref_img->height - r - 1; y += grid_step)
			for (int x = r + 1; x < ref_img->width - r - 1; x += grid_step)
			{
				if (on_ground && on_ground->validCount(x, y, r, r) < side * side) continue;
				POI2D poi(x, y);
				if (perspective) coarse_h.seed(poi);
				else coarse.seed(poi);
				const Point2D t(x + poi.deformation.u, y + poi.deformation.v);
				if (t.x < r + 1 || t.y < r + 1 || t.x > tar_img->width - r - 2 || t.y > tar_img->height - r - 2) continue;
				queue.push_back(poi);
			}
		result.refined_points = (int)queue.size();

		result.status = AlignmentStatus::too_few_ground_points;
		if (result.refined_points < min_refined) return result;

		//⚑ Under a tilt the camera's own warp across one subset is not affine, and a
		//first-order subset finds its centre a little off where the warp curves, in a
		//pattern that varies over the picture and so moves the fit: 0.021 px worst on the
		//smoke test's 5 % tilt, against 0.003 px at the subsets themselves. A second-order
		//subset follows the curvature, so the perspective refinement uses one.
		if (perspective)
		{
			ICGN2D2 icgn2(r, r, 0.001f, 20, thread_number);
			icgn2.setImages(*ref_img, *tar_img);
			icgn2.prepare();
			icgn2.compute(queue);
		}
		else
			icgn.compute(queue);

		std::vector<Point2D> centre, found;
		std::vector<size_t> origin; //centre[k] is queue[origin[k]]
		for (size_t i = 0; i < queue.size(); i++)
		{
			const POI2D& poi = queue[i];
			result.refined_centre.push_back(Point2D(poi.x, poi.y));
			result.refined_zncc.push_back(poi.result.zncc);
			result.refined_residual.push_back(std::numeric_limits<float>::quiet_NaN());
			result.refined_kept.push_back(0);
			if (!(poi.result.zncc >= 0.9f) || !std::isfinite(poi.deformation.u) || !std::isfinite(poi.deformation.v)) continue;
			centre.push_back(Point2D(poi.x, poi.y));
			found.push_back(Point2D(poi.x + poi.deformation.u, poi.y + poi.deformation.v));
			origin.push_back(i);
		}
		if ((int)centre.size() < min_refined) return result;

		//the similarity, and the homography through the same subsets: one is the camera
		//move, the other a check on it (or, asked for perspective, the other way round)
		Similarity2D fine = coarse;
		std::vector<char> keep(centre.size(), 1);
		if (fixed_scale > 0 && !fitSimilarity(centre, found, keep, fine, fixed_scale)) return result;
		const bool similar = reweighted(centre, found, keep, fine,
			[&](const std::vector<char>& k, Similarity2D& out) { return fitSimilarity(centre, found, k, out, fixed_scale); });
		if (!similar && !perspective) return result;

		Homography2D projective = perspective ? coarse_h : Homography2D::from(fine);
		std::vector<char> keep_h(centre.size(), 1);
		result.projective_fitted = fitHomography(centre, found, keep_h, projective)
			&& reweighted(centre, found, keep_h, projective,
				[&](const std::vector<char>& k, Homography2D& out) { return fitHomography(centre, found, k, out); });
		if (perspective && !result.projective_fitted) return result;

		auto account = [&](const auto& camera, const std::vector<char>& kept, bool record, int& used) {
			double ss = 0;
			used = 0;
			for (size_t i = 0; i < centre.size(); i++)
			{
				const double e = residual(camera, centre[i], found[i]);
				if (record)
				{
					result.refined_residual[origin[i]] = (float)e;
					result.refined_kept[origin[i]] = kept[i];
				}
				if (!kept[i]) continue;
				ss += e * e;
				used++;
			}
			return used > 0 ? (float)std::sqrt(ss / used) : 0.f;
		};
		int used = 0, used_h = 0;
		const float rms = similar ? account(fine, keep, !perspective, used) : 0.f;
		if (result.projective_fitted)
		{
			result.projective = projective;
			result.projective_rms = account(projective, keep_h, perspective, used_h);
		}
		if (similar) result.transform = fine;

		result.perspective = perspective;
		result.refined_used = perspective ? used_h : used;
		if (result.refined_used < min_refined) return result;
		result.residual_rms = perspective ? result.projective_rms : rms;
		result.status = AlignmentStatus::ok;
		return result;
	}

}//namespace opencorr
