/*
 Local build-verification test for SessionAlignment2D: a reference and a target
 photographed in separate sessions, so the camera has moved between them and
 the target need not even be the same size.

 The scene is rendered analytically, so the answer is exact. A static rig is
 covered in speckle, with a specimen disc on it that deforms between sessions
 by d(X) = d0 + E (X - c). The reference camera sees world coordinates
 directly; the target camera maps a world point W to x = s R W + t, in an image
 of a different size. Every coefficient is lopsided -- E has four different
 non-zero terms, the camera both scales and rotates, the shift differs in x and
 y -- so no two quantities the code could confuse happen to agree.

 What is held, and why each is a separate case:

 1. Fitted on the GROUND (the rig, marked as a region), the camera move is
    recovered to within a few thousandths of a pixel anywhere on the image.
    The tolerance is set below what SIFT and RANSAC reach alone, so the
    sub-pixel refinement is what this case pins.
 2. A specimen point measured with the alignment as its starting guess, and
    carried back into the reference frame, reports the specimen's own
    deformation: d0 and all four terms of E.
 3. With no ground marked, the whole image is used. The specimen then makes up
    a third of the matches and moves differently from the rest, and the
    refinement must reject it rather than average it in.
 4. A target with nothing to match is refused with its own status, never
    answered with an identity transform that reads as "the camera did not
    move".
 5. Similarity2D's own arithmetic, on a point whose answer is worked by hand.
 6. A specimen that is MOST of the picture. Rejection alone then keeps the
    specimen's motion as the consensus, so with no ground marked the fit is
    1.3 px wrong, and a narrow marked strip gets it right. This is the case that
    says what the ground is for; in the scene above, the rejection removes the
    specimen either way.
 7. The same, with the specimen slid 9 px: beyond what the refinement reaches
    from a bad starting fit, so the ground must already choose the SIFT matches.
 8. No ground, fitted on the specimen, with a KNOWN LENGTH between two rig
    points: the scale is taken from it, so the specimen's uniform dilatation is
    measured. The same fit without it loses the dilatation from both normal
    terms, which is what lets the case tell the two apart.
 9. A known length with an end too near the edge is refused.
10. A known length with an end covered on the target is refused.

 Measured: ground fit 0.0012 px worst over the image (RANSAC alone 0.033 px);
 specimen displacement carried back within 0.0011 px on average.

 Negative checks, each run and reverted:
 - Reporting the RANSAC fit instead of the refined one: cases 1, 2, 3 red
   (0.033 px, and 0.052 px with no ground).
 - No rejection in the refinement: case 3 red, 0.64 px.
 - toReference without the camera's rotation: cases 2 and 5 red.
 - unmap with the rotation transposed: cases 2 and 5 red.
 - seed with the sign of uy flipped: case 5 red only. ICGN converges from the
   wrong gradient anyway on this pattern, so the alignment cases cannot see it.
 - The refinement grid ignoring the ground: case 6 red, 1.3 px.
 - SIFT ignoring the ground: case 7 red (refused, no ground point converged).
   Cases 1 to 6 stay green: the refinement recovers from a coarse fit 1.3 px off.
 - A ground subset admitted when only its CENTRE is on the ground: NOT caught.
   The rejection removes the subsets that straddle onto the specimen (267 of 321
   kept). The whole-subset rule stays, so the fit does not depend on that.
 - SIFT2D's empty-descriptor guard removed: case 4 aborts in FLANN.

 The known length (2026-10-16), each run and reverted:
 - The fixed scale ignored by the fit, or the ratio of the lengths inverted, or
   the final rounds fitted freely: case 8 red on the scale and the strains.
 - The correlation of an end not checked: case 10 red.
 - The bounds of an end not checked: NOT caught. ICGN refuses a subset outside
   the image itself, so case 9 is refused either way; the check stays so that
   nothing is read past an image edge to find out.
 - The first fit at the known scale, before the reweighting, removed: NOT
   caught. The rounds that follow reach the same fit from the free one.
*/

#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

namespace
{
	int failures = 0;

	void check(bool ok, const char* what)
	{
		printf("%s  %s\n", ok ? "PASS" : "FAIL", what);
		if (!ok) failures++;
	}

	struct Speckle
	{
		int cells_x, cells_y;
		float x0, y0;
		vector<float> dot_x, dot_y;
		static constexpr float cell = 6.f, radius = 2.2f;

		Speckle(float x0, float y0, int width, int height, unsigned seed) : x0(x0), y0(y0)
		{
			cells_x = (int)(width / cell) + 8;
			cells_y = (int)(height / cell) + 8;
			mt19937 rng(seed);
			uniform_real_distribution<float> jitter(0.f, cell);
			for (int j = 0; j < cells_y; j++)
				for (int i = 0; i < cells_x; i++)
				{
					dot_x.push_back(x0 + (i - 4) * cell + jitter(rng));
					dot_y.push_back(y0 + (j - 4) * cell + jitter(rng));
				}
		}

		float at(float x, float y) const
		{
			const int ci = (int)floor((x - x0) / cell) + 4;
			const int cj = (int)floor((y - y0) / cell) + 4;
			float sum = 0.f;
			for (int j = cj - 3; j <= cj + 3; j++)
				for (int i = ci - 3; i <= ci + 3; i++)
				{
					if (i < 0 || j < 0 || i >= cells_x || j >= cells_y) continue;
					const float dx = x - dot_x[j * cells_x + i], dy = y - dot_y[j * cells_x + i];
					sum += exp(-(dx * dx + dy * dy) / (radius * radius));
				}
			return 20.f + 180.f * min(sum, 1.f);
		}
	};

	const float CX = 120.f, CY = 100.f;
	float DISC = 60.f; //case 6 enlarges it until the specimen is most of the picture
	const double E11 = 0.012, E12 = 0.004, E21 = 0.002, E22 = -0.005;
	double D0X = 1.3, D0Y = -0.7; //case 7 slides the specimen much further
	const double SCALE = 1.07, ANGLE = 3.5 * M_PI / 180.0;

	Point2D deform(double X, double Y) //material -> world after deformation
	{
		return Point2D((float)(X + D0X + E11 * (X - CX) + E12 * (Y - CY)),
			(float)(Y + D0Y + E21 * (X - CX) + E22 * (Y - CY)));
	}

	Point2D undeform(double x, double y) //world after deformation -> material
	{
		const double a = 1 + E11, b = E12, c = E21, d = 1 + E22, det = a * d - b * c;
		const double rx = x - CX - D0X, ry = y - CY - D0Y;
		return Point2D((float)(CX + (d * rx - b * ry) / det), (float)(CY + (-c * rx + a * ry) / det));
	}

	//the true camera move, with its own arithmetic: built from the class under test, a
	//broken map() would render the target and judge the answer alike, and agree with
	//itself (it did, against a stub returning the identity)
	struct Camera
	{
		double s, theta, tx, ty;
		Point2D map(double X, double Y) const
		{
			return Point2D((float)(s * (cos(theta) * X - sin(theta) * Y) + tx),
				(float)(s * (sin(theta) * X + cos(theta) * Y) + ty));
		}
		Point2D unmap(double x, double y) const
		{
			const double rx = (x - tx) / s, ry = (y - ty) / s;
			return Point2D((float)(cos(theta) * rx + sin(theta) * ry), (float)(-sin(theta) * rx + cos(theta) * ry));
		}
	};

	bool inDisc(double X, double Y) { return (X - CX) * (X - CX) + (Y - CY) * (Y - CY) <= DISC * DISC; }

	//camera is null for the reference, which sees the world directly
	Image2D render(int w, int h, const Camera* camera, bool deformed, const Speckle& specimen, const Speckle& rig)
	{
		Image2D img(w, h);
		const int ss = 3;
		for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++)
			{
				float acc = 0.f;
				for (int sy = 0; sy < ss; sy++)
					for (int sx = 0; sx < ss; sx++)
					{
						const Point2D p(x + (sx + 0.5f) / ss - 0.5f, y + (sy + 0.5f) / ss - 0.5f);
						const Point2D W = camera ? camera->unmap(p.x, p.y) : p;
						const Point2D X = deformed ? undeform(W.x, W.y) : W;
						acc += inDisc(X.x, X.y) ? specimen.at(X.x, X.y) : rig.at(W.x, W.y);
					}
				img.eg_mat(y, x) = acc / (ss * ss);
				img.cv_mat.at<uchar>(y, x) = (uchar)lround(min(255.f, acc / (ss * ss)));
			}
		return img;
	}

	//the largest distance between where the fitted and the true transform put a
	//reference pixel, over the whole reference image
	double worstMapping(const Similarity2D& fitted, const Camera& truth, int w, int h)
	{
		double worst = 0;
		for (int y = 0; y < h; y += 4)
			for (int x = 0; x < w; x += 4)
			{
				const Point2D a = fitted.map(Point2D((float)x, (float)y));
				const Point2D b = truth.map(x, y);
				worst = max(worst, (double)hypot(a.x - b.x, a.y - b.y));
			}
		return worst;
	}

	unique_ptr<Shape2D> rectangle(int x0, int y0, int x1, int y1)
	{
		return make_unique<Polygon2D>(vector<int>{ x0, x1, x1, x0 }, vector<int>{ y0, y0, y1, y1 });
	}
}

int main()
{
	const int ref_w = 240, ref_h = 200, tar_w = 280, tar_h = 240;
	const Speckle specimen(-20, -20, 300, 260, 11), rig(-80, -80, 420, 380, 23);

	//the target camera, placed so the specimen centre lands at (140, 120)
	Camera truth{ SCALE, ANGLE, 0, 0 };
	{
		const Point2D c = truth.map(CX, CY);
		truth.tx = 140 - c.x;
		truth.ty = 120 - c.y;
	}

	Image2D ref = render(ref_w, ref_h, nullptr, false, specimen, rig);
	Image2D tar = render(tar_w, tar_h, &truth, true, specimen, rig);

	//the ground: the whole reference less the disc and a 10 px margin around it
	OrderedRegion2D ground;
	ground.add(rectangle(0, 0, ref_w - 1, ref_h - 1));
	ground.cut(make_unique<Ellipse2D>(CX, CY, DISC + 10, DISC + 10));

	// 1. fitted on the ground
	SessionAlignment2D alignment(15, 10, 4);
	alignment.setImages(ref, tar);
	alignment.setGround(&ground);
	const AlignmentResult on_ground = alignment.compute();
	{
		const double worst = worstMapping(on_ground.transform, truth, ref_w, ref_h);
		printf("ground: %d matches, %d inliers, %d of %d refined points used, rms %.4f px, worst mapping error %.5f px\n",
			on_ground.matches, on_ground.inliers, on_ground.refined_used, on_ground.refined_points,
			on_ground.residual_rms, worst);
		printf("        scale %.6f (true %.6f), angle %.5f deg (true %.5f)\n", on_ground.transform.scale, SCALE,
			on_ground.transform.angle * 180 / M_PI, ANGLE * 180 / M_PI);
		check(on_ground.status == AlignmentStatus::ok, "the ground alignment succeeds");
		check(worst < 0.005, "the ground alignment places every reference pixel within 0.005 px of the truth");
		check(on_ground.refined_used >= 20, "the refinement rests on at least 20 ground subsets");
		check(on_ground.residual_rms < 0.01f, "the refined ground points sit within 0.01 px rms of the fit");
	}

	// 2. a specimen point measured from the alignment and carried back
	{
		ICGN2D1 icgn(16, 16, 0.001f, 20, 4);
		icgn.setImages(ref, tar);
		icgn.prepare();
		vector<POI2D> pois;
		for (int y = (int)CY - 36; y <= (int)CY + 36; y += 6)
			for (int x = (int)CX - 36; x <= (int)CX + 36; x += 6)
			{
				if ((x - CX) * (x - CX) + (y - CY) * (y - CY) > 40 * 40) continue;
				POI2D poi(Point2D((float)x, (float)y));
				on_ground.transform.seed(poi);
				pois.push_back(poi);
			}
		icgn.compute(pois);

		//a single subset's first-order gradient scatters by about 1e-3 on this pattern,
		//so the case holds each term's MEAN: a wrong conversion is off by the camera's own
		//term, about 0.07, in every point alike
		double error = 0, mean_gradient[4] = { 0, 0, 0, 0 };
		int solved = 0;
		for (auto& poi : pois)
		{
			if (poi.result.zncc < 0.9f) continue;
			on_ground.transform.toReference(poi);
			const Point2D t = deform(poi.x, poi.y);
			error += hypot(poi.deformation.u - (t.x - poi.x), poi.deformation.v - (t.y - poi.y));
			mean_gradient[0] += poi.deformation.ux;
			mean_gradient[1] += poi.deformation.uy;
			mean_gradient[2] += poi.deformation.vx;
			mean_gradient[3] += poi.deformation.vy;
			solved++;
		}
		error /= max(solved, 1);
		const double truth_gradient[4] = { E11, E12, E21, E22 };
		double worst_gradient = 0;
		for (int k = 0; k < 4; k++)
		{
			mean_gradient[k] /= max(solved, 1);
			worst_gradient = max(worst_gradient, fabs(mean_gradient[k] - truth_gradient[k]));
		}
		printf("specimen: %d of %zu solved, mean displacement error %.5f px\n", solved, pois.size(), error);
		printf("          mean gradient %+.5f %+.5f %+.5f %+.5f (true %+.5f %+.5f %+.5f %+.5f)\n",
			mean_gradient[0], mean_gradient[1], mean_gradient[2], mean_gradient[3], E11, E12, E21, E22);
		check(solved == (int)pois.size(), "every specimen point solves from the alignment's starting guess");
		check(error < 0.01, "the specimen's displacement, carried back, is within 0.01 px of the truth");
		check(solved > 0 && worst_gradient < 2e-4, "each displacement gradient term, carried back, averages within 2e-4 of E");
	}

	// 3. no ground marked: the specimen must be rejected, not averaged in
	{
		SessionAlignment2D whole(15, 10, 4);
		whole.setImages(ref, tar);
		const AlignmentResult r = whole.compute();
		const double worst = worstMapping(r.transform, truth, ref_w, ref_h);
		printf("whole image: %d matches, %d inliers, %d of %d refined points used, worst mapping error %.5f px\n",
			r.matches, r.inliers, r.refined_used, r.refined_points, worst);
		check(r.status == AlignmentStatus::ok, "the whole-image alignment succeeds");
		check(worst < 0.005, "with no ground marked the moving specimen is rejected, not averaged in");
		check(r.refined_used < r.refined_points, "the whole-image refinement rejected some points");
	}

	// 4. nothing to match
	{
		Image2D blank(tar_w, tar_h);
		blank.cv_mat.setTo(128);
		blank.eg_mat.setConstant(128.f);
		SessionAlignment2D none(15, 10, 4);
		none.setImages(ref, blank);
		const AlignmentResult r = none.compute();
		check(r.status == AlignmentStatus::too_few_matches, "a featureless target is refused as too few matches");
	}

	// 5. Similarity2D arithmetic, by hand: scale 2, a quarter turn, shift (10, 20)
	{
		Similarity2D s;
		s.scale = 2;
		s.angle = M_PI / 2;
		s.tx = 10;
		s.ty = 20;
		const Point2D m = s.map(Point2D(3, 1)); //2 * (-1, 3) + (10, 20) = (8, 26)
		check(fabs(m.x - 8) < 1e-5 && fabs(m.y - 26) < 1e-5, "map: (3, 1) goes to (8, 26)");
		const Point2D u = s.unmap(m);
		check(fabs(u.x - 3) < 1e-5 && fabs(u.y - 1) < 1e-5, "unmap undoes map");

		POI2D poi(Point2D(3, 1));
		s.seed(poi);
		check(fabs(poi.deformation.u - 5) < 1e-5 && fabs(poi.deformation.v - 25) < 1e-5,
			"seed: the starting displacement is where map puts the point");
		check(fabs(poi.deformation.ux + 1) < 1e-5 && fabs(poi.deformation.uy + 2) < 1e-5
			&& fabs(poi.deformation.vx - 2) < 1e-5 && fabs(poi.deformation.vy + 1) < 1e-5,
			"seed: the starting gradient is the camera's own (2R - I)");

		//a target-frame measurement of a point that in fact moved by (0.5, -0.25) and was
		//stretched by 1 % in x: carried back, exactly that is reported
		POI2D measured(Point2D(3, 1));
		//(3.5, 0.75) turned a quarter is (-0.75, 3.5), doubled (-1.5, 7), shifted (8.5, 27)
		measured.deformation.u = 8.5f - 3;
		measured.deformation.v = 27.f - 1;
		//G = s R F, F = [[1.01, 0], [0, 1]]: s R = [[0, -2], [2, 0]]
		measured.deformation.ux = 0 - 1;
		measured.deformation.uy = -2;
		measured.deformation.vx = 2 * 1.01f;
		measured.deformation.vy = 0 - 1;
		s.toReference(measured);
		check(fabs(measured.deformation.u - 0.5) < 1e-5 && fabs(measured.deformation.v + 0.25) < 1e-5,
			"toReference: the displacement comes back in reference pixels");
		check(fabs(measured.deformation.ux - 0.01) < 1e-5 && fabs(measured.deformation.uy) < 1e-5
			&& fabs(measured.deformation.vx) < 1e-5 && fabs(measured.deformation.vy) < 1e-5,
			"toReference: the gradient comes back as the specimen's own");
	}

	// 6. a specimen that is MOST of the picture: the rejection alone now keeps the
	//    specimen's motion as the consensus and throws the rig away, and only the
	//    marked ground gets the camera right
	{
		DISC = 92.f;
		Image2D big_ref = render(ref_w, ref_h, nullptr, false, specimen, rig);
		Image2D big_tar = render(tar_w, tar_h, &truth, true, specimen, rig);
		OrderedRegion2D strip;
		strip.add(rectangle(0, 0, ref_w - 1, ref_h - 1));
		strip.cut(make_unique<Ellipse2D>(CX, CY, DISC + 6, DISC + 6));

		SessionAlignment2D with_ground(8, 6, 4);
		with_ground.setImages(big_ref, big_tar);
		with_ground.setGround(&strip);
		const AlignmentResult g = with_ground.compute();

		SessionAlignment2D without(8, 6, 4);
		without.setImages(big_ref, big_tar);
		const AlignmentResult w = without.compute();

		const double worst_g = worstMapping(g.transform, truth, ref_w, ref_h);
		const double worst_w = worstMapping(w.transform, truth, ref_w, ref_h);
		printf("mostly specimen, ground marked: %d matches, %d of %d refined used, worst mapping error %.5f px\n",
			g.matches, g.refined_used, g.refined_points, worst_g);
		printf("mostly specimen, none marked:   %d matches, %d of %d refined used, worst mapping error %.5f px\n",
			w.matches, w.refined_used, w.refined_points, worst_w);
		check(g.status == AlignmentStatus::ok && worst_g < 0.01, "a narrow strip of marked ground still gets the camera right");
		check(worst_w > 0.1, "with no ground marked, a specimen that is most of the picture is taken for the camera");
		DISC = 60.f;
	}

	// 7. the same, with the specimen slid 9 px in its grips: the RANSAC fit on
	//    everything is now 9 px wrong on the ground, beyond what the refinement can
	//    reach from, so the ground must already govern which matches are fitted
	{
		DISC = 92.f;
		D0X = 7.0;
		D0Y = -5.5;
		Image2D big_ref = render(ref_w, ref_h, nullptr, false, specimen, rig);
		Image2D big_tar = render(tar_w, tar_h, &truth, true, specimen, rig);
		OrderedRegion2D strip;
		strip.add(rectangle(0, 0, ref_w - 1, ref_h - 1));
		strip.cut(make_unique<Ellipse2D>(CX, CY, DISC + 6, DISC + 6));

		SessionAlignment2D with_ground(8, 6, 4);
		with_ground.setImages(big_ref, big_tar);
		with_ground.setGround(&strip);
		const AlignmentResult g = with_ground.compute();
		const double worst_g = worstMapping(g.transform, truth, ref_w, ref_h);
		printf("slid specimen, ground marked: %d matches, %d inliers, %d of %d refined used, worst mapping error %.5f px\n",
			g.matches, g.inliers, g.refined_used, g.refined_points, worst_g);
		check(g.status == AlignmentStatus::ok && worst_g < 0.01, "a specimen slid 9 px does not pull the ground fit off");
		DISC = 60.f;
		D0X = 1.3;
		D0Y = -0.7;
	}

	// 8. no ground, fitted on the specimen, with a KNOWN LENGTH: two points on the rig
	//    whose distance did not change. The camera's scale is taken from them, so the
	//    fit on the specimen absorbs only its rigid motion and its uniform dilatation is
	//    measured. Without the length the same fit loses (E11 + E22) / 2 from both
	//    normal terms, which is what makes this case able to tell the two apart.
	{
		OrderedRegion2D disc;
		disc.add(make_unique<Ellipse2D>(CX, CY, DISC - 12, DISC - 12));

		auto specimenGradient = [&](const AlignmentResult& r, double out[4]) {
			ICGN2D1 icgn(16, 16, 0.001f, 20, 4);
			icgn.setImages(ref, tar);
			icgn.prepare();
			vector<POI2D> pois;
			for (int y = (int)CY - 36; y <= (int)CY + 36; y += 6)
				for (int x = (int)CX - 36; x <= (int)CX + 36; x += 6)
				{
					if ((x - CX) * (x - CX) + (y - CY) * (y - CY) > 40 * 40) continue;
					POI2D poi(Point2D((float)x, (float)y));
					r.transform.seed(poi);
					pois.push_back(poi);
				}
			icgn.compute(pois);
			int solved = 0;
			for (int k = 0; k < 4; k++) out[k] = 0;
			for (auto& poi : pois)
			{
				if (poi.result.zncc < 0.9f) continue;
				r.transform.toReference(poi);
				out[0] += poi.deformation.ux;
				out[1] += poi.deformation.uy;
				out[2] += poi.deformation.vx;
				out[3] += poi.deformation.vy;
				solved++;
			}
			for (int k = 0; k < 4; k++) out[k] /= max(solved, 1);
			return solved;
		};

		SessionAlignment2D plain(15, 6, 4);
		plain.setImages(ref, tar);
		plain.setGround(&disc);
		const AlignmentResult p = plain.compute();

		SessionAlignment2D known(15, 6, 4);
		known.setImages(ref, tar);
		known.setGround(&disc);
		known.setKnownLength(Point2D(22, 24), Point2D(214, 178));
		const AlignmentResult k = known.compute();

		double gp[4], gk[4];
		const int np = specimenGradient(p, gp), nk = specimenGradient(k, gk);
		printf("specimen fit, no length:   scale %.6f, %d points, gradient %+.5f %+.5f %+.5f %+.5f\n",
			p.transform.scale, np, gp[0], gp[1], gp[2], gp[3]);
		printf("specimen fit, known length: scale %.6f, %d points, gradient %+.5f %+.5f %+.5f %+.5f\n",
			k.transform.scale, nk, gk[0], gk[1], gk[2], gk[3]);
		printf("                            ends found at (%.3f, %.3f) and (%.3f, %.3f), zncc %.4f %.4f\n",
			k.known_length_target[0].x, k.known_length_target[0].y, k.known_length_target[1].x,
			k.known_length_target[1].y, k.known_length_zncc[0], k.known_length_zncc[1]);
		const Point2D ta = truth.map(22, 24), tb = truth.map(214, 178);
		check(k.status == AlignmentStatus::ok && k.known_length_used, "the known length is measured and used");
		check(hypot(k.known_length_target[0].x - ta.x, k.known_length_target[0].y - ta.y) < 0.01
			&& hypot(k.known_length_target[1].x - tb.x, k.known_length_target[1].y - tb.y) < 0.01,
			"each end of the known length is found on the target within 0.01 px");
		check(fabs(k.transform.scale - SCALE) < 1e-4, "the camera scale is taken from the known length");
		check(!p.known_length_used && fabs(gp[0] - E11) > 2e-3 && fabs(gp[3] - E22) > 2e-3,
			"without the length, the specimen fit loses the dilatation from both normal terms");
		//the rotation is still the fit's: only the symmetric part is the specimen's
		check(nk > 0 && fabs(gk[0] - E11) < 3e-4 && fabs(gk[3] - E22) < 3e-4
			&& fabs((gk[1] + gk[2]) - (E12 + E21)) < 3e-4,
			"with it, both normal strains and the shear are measured in full");
	}

	// 9. a known length whose end cannot be measured is refused, never quietly ignored
	{
		SessionAlignment2D edge(15, 10, 4);
		edge.setImages(ref, tar);
		edge.setGround(&ground);
		edge.setKnownLength(Point2D(3, 3), Point2D(214, 178));
		const AlignmentResult r = edge.compute();
		check(r.status == AlignmentStatus::known_length_unmeasured && !r.known_length_used,
			"a known length with an end too near the edge to correlate is refused");
	}

	// 10. an end inside both images that does not correlate: the target is covered
	//     where it lands, as by a hand or a cable, so only the correlation can refuse it
	{
		Image2D covered = render(tar_w, tar_h, &truth, true, specimen, rig);
		const Point2D tb = truth.map(214, 178);
		for (int y = (int)tb.y - 20; y <= (int)tb.y + 20; y++)
			for (int x = (int)tb.x - 20; x <= (int)tb.x + 20; x++)
			{
				if (x < 0 || y < 0 || x >= tar_w || y >= tar_h) continue;
				covered.eg_mat(y, x) = 128.f + 40.f * (float)((x / 7 + y / 5) % 2);
				covered.cv_mat.at<uchar>(y, x) = (uchar)covered.eg_mat(y, x);
			}
		SessionAlignment2D hidden(15, 10, 4);
		hidden.setImages(ref, covered);
		hidden.setGround(&ground);
		hidden.setKnownLength(Point2D(22, 24), Point2D(214, 178));
		const AlignmentResult r = hidden.compute();
		printf("covered end: zncc %.4f %.4f\n", r.known_length_zncc[0], r.known_length_zncc[1]);
		check(r.status == AlignmentStatus::known_length_unmeasured && !r.known_length_used,
			"a known length with an end covered on the target is refused");
	}

	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}
