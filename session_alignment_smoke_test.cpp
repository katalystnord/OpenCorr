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
11. Ground marked partly OVER THE SPECIMEN, with a patch of real ground covered
    on the target. Every refinement subset is reported with its centre, its
    correlation, its offset from the fit and whether the fit kept it, so a
    caller can say where ground moved. The specimen's subsets are rejected 1.26
    px or more from the fit, the kept ones sit within 0.003 px.
12. Homography2D's arithmetic by hand: map, unmap, the Jacobian, and a point
    that moved only with the camera carrying back to no movement.
13. The camera TILTED between sessions (perspective terms giving about 5 %
    scale change across the picture). Asked for perspective, the fit places
    every reference pixel within 0.002 px and the specimen reads its own E;
    the similarity cannot (16.9 px). Unasked, the perspective CHECK fitted
    through the same subsets departs from the similarity by pixels, so a caller
    can tell a tilt from one run.
14. On the untilted pair the check agrees with the similarity within 0.0014 px.

 ⚑ Under a tilt the camera's own warp across one subset curves. First-order
 subsets find their centres a little off in a pattern that varies over the
 picture: the perspective fit was 0.021 px worst with them, 0.0017 px with
 second-order ones, which the refinement now uses. The specimen MEASURED with
 first-order subsets through a correct fit is 0.023 px off on average for the
 same reason, 0.0023 px with second-order ones; that is the caller's to say.

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

 Every refinement subset reported (2026-10-17), each run and reverted:
 - Offsets and kept flags indexed by the CONVERGED subsets instead of all of
   them: case 11 red, four checks. ⚑ Not caught until the covered patch was
   added: with every subset converging the two indexings are the same.
 - An offset written only for kept subsets: case 11 red (also not caught until
   the offsets were required to be finite: a min() over NaN keeps its start).
 - Unconverged subsets left out of the report: case 11 red, six checks.
 - Every subset reported kept: case 11 red, three checks.
 - The offset measured from the RANSAC fit, not the final one: case 11 red,
   and case 1.
 - A subset that did not correlate given an offset of zero: case 11 red.

 Perspective (2026-10-17), each run and reverted:
 - toReference taking the Jacobian at the reference point, not where the
   material point now sits: case 13 red on the gradient.
 - toReference not inverting the Jacobian: cases 12 and 13 red.
 - The Jacobian without its perspective terms, or with a sign flipped: 12, 13.
 - unmap with the adjugate transposed: 12 and 13 red, four checks.
 - The refinement with first-order subsets: case 13 red (0.021 px).
 - The perspective flag not reported, or the check not fitted unasked: 13, 14.
 - The residual reported from the similarity when perspective was asked: 13.
 - The refinement seeded from the similarity instead of the RANSAC homography:
   case 13 red, subsets at the edges no longer converge (counts fall short).
 - The homography fit not reweighted: case 11 red, the check pulled 1.3 px
   toward the ground that moved.
 - Gauss-Newton cut to one iteration: NOT caught. In normalised coordinates
   the fit is nearly linear and one step from a similarity start lands within
   the tolerances here; the iterations stay for stronger tilts.
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

	//a camera TILTED against the reference: x = H W, with perspective terms, so the
	//target's scale varies across the picture. Its own arithmetic, as Camera has.
	struct TiltedCamera
	{
		double h[9];
		Point2D map(double X, double Y) const
		{
			const double w = h[6] * X + h[7] * Y + h[8];
			return Point2D((float)((h[0] * X + h[1] * Y + h[2]) / w), (float)((h[3] * X + h[4] * Y + h[5]) / w));
		}
		Point2D unmap(double x, double y) const
		{
			//the adjugate of h: the inverse up to a scale the division removes
			const double a = h[4] * h[8] - h[5] * h[7], b = h[2] * h[7] - h[1] * h[8], c = h[1] * h[5] - h[2] * h[4];
			const double d = h[5] * h[6] - h[3] * h[8], e = h[0] * h[8] - h[2] * h[6], f = h[2] * h[3] - h[0] * h[5];
			const double g = h[3] * h[7] - h[4] * h[6], k = h[1] * h[6] - h[0] * h[7], l = h[0] * h[4] - h[1] * h[3];
			const double w = g * x + k * y + l;
			return Point2D((float)((a * x + b * y + c) / w), (float)((d * x + e * y + f) / w));
		}
	};

	bool inDisc(double X, double Y) { return (X - CX) * (X - CX) + (Y - CY) * (Y - CY) <= DISC * DISC; }

	//camera is null for the reference, which sees the world directly
	Image2D render(int w, int h, const Camera* camera, bool deformed, const Speckle& specimen, const Speckle& rig,
		const TiltedCamera* tilted = nullptr)
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
						const Point2D W = tilted ? tilted->unmap(p.x, p.y) : camera ? camera->unmap(p.x, p.y) : p;
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

	template <class A, class B>
	double worstAgainst(const A& a, const B& b, int w, int h)
	{
		double worst = 0;
		for (int y = 0; y < h; y += 4)
			for (int x = 0; x < w; x += 4)
			{
				const Point2D p = a.map(Point2D((float)x, (float)y)), q = b.map(Point2D((float)x, (float)y));
				worst = max(worst, (double)hypot(p.x - q.x, p.y - q.y));
			}
		return worst;
	}

	template <class Fitted>
	double worstAgainstTilted(const Fitted& fitted, const TiltedCamera& truth, int w, int h)
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

	// 11. ground marked partly OVER THE SPECIMEN: the rig, plus a box on the disc that
	//     the caller wrongly calls stationary. The refinement rejects the box's
	//     subsets, and the result must say which subsets those were, where they sit
	//     and how far they disagree with the fit, so the caller can say WHERE the
	//     ground moved rather than only that a count fell short.
	{
		OrderedRegion2D wrong;
		wrong.add(rectangle(0, 0, ref_w - 1, ref_h - 1));
		wrong.cut(make_unique<Ellipse2D>(CX, CY, DISC + 10, DISC + 10));
		const int bx0 = (int)CX - 40, by0 = (int)CY - 30, bx1 = (int)CX + 10, by1 = (int)CY + 30;
		wrong.add(rectangle(bx0, by0, bx1, by1)); //inside the disc: this moved

		//and a patch of real ground covered on the target, as by a cable, so some
		//subsets do not correlate at all: without them every subset converges, and a
		//report indexed by the converged subsets lines up with one indexed by all of them
		const int cx0 = 175, cy0 = 20, cx1 = 230, cy1 = 60; //reference pixels, on the rig
		Image2D covered = render(tar_w, tar_h, &truth, true, specimen, rig);
		for (int y = 0; y < tar_h; y++)
			for (int x = 0; x < tar_w; x++)
			{
				const Point2D W = truth.unmap(x, y);
				if (W.x < cx0 || W.x > cx1 || W.y < cy0 || W.y > cy1) continue;
				covered.eg_mat(y, x) = 90.f;
				covered.cv_mat.at<uchar>(y, x) = 90;
			}

		SessionAlignment2D partly(15, 10, 4);
		partly.setImages(ref, covered);
		partly.setGround(&wrong);
		const AlignmentResult r = partly.compute();
		const double worst = worstMapping(r.transform, truth, ref_w, ref_h);

		const size_t n = r.refined_centre.size();
		const bool sized = n == (size_t)r.refined_points && r.refined_zncc.size() == n
			&& r.refined_residual.size() == n && r.refined_kept.size() == n;
		int kept = 0, in_box = 0, in_box_rejected = 0, rejected_off_box = 0, unmatched = 0, unmatched_reported = 0;
		float least_box_residual = 1e9f, worst_kept_residual = 0.f;
		bool box_offsets_finite = true;
		if (sized)
			for (size_t i = 0; i < n; i++)
			{
				const Point2D& c = r.refined_centre[i];
				const bool box = c.x >= bx0 && c.x <= bx1 && c.y >= by0 && c.y <= by1;
				if (r.refined_kept[i])
				{
					kept++;
					worst_kept_residual = max(worst_kept_residual, r.refined_residual[i]);
				}
				if (!(r.refined_zncc[i] >= 0.9f))
				{
					unmatched++;
					if (!r.refined_kept[i] && std::isnan(r.refined_residual[i])) unmatched_reported++;
					continue;
				}
				if (box)
				{
					in_box++;
					if (!r.refined_kept[i]) in_box_rejected++;
					if (!std::isfinite(r.refined_residual[i])) box_offsets_finite = false;
					else least_box_residual = min(least_box_residual, r.refined_residual[i]);
				}
				else if (!r.refined_kept[i] && r.refined_residual[i] > 0.1f)
					rejected_off_box++;
			}
		printf("ground partly on the specimen: %d of %d used, %d box subsets of which %d rejected, "
			"least box offset %.3f px, worst kept %.4f px, %d not correlated, worst mapping error %.5f px\n",
			r.refined_used, r.refined_points, in_box, in_box_rejected, least_box_residual,
			worst_kept_residual, unmatched, worst);
		check(r.status == AlignmentStatus::ok && worst < 0.005, "ground partly on the specimen still gets the camera right");
		check(sized, "every refinement subset is reported, with its correlation, its offset and whether it was kept");
		check(sized && kept == r.refined_used, "the kept subsets reported are the ones the fit used");
		check(in_box > 0 && in_box_rejected == in_box, "every ground subset on the moving specimen is reported rejected");
		check(in_box > 0 && box_offsets_finite && least_box_residual > 0.5f,
			"and each is reported as far from the fit as the specimen moved");
		check(unmatched > 0 && unmatched_reported == unmatched,
			"a subset that did not correlate is reported unkept, with no offset");
		check(rejected_off_box == 0, "no subset on the real ground is reported as having moved by more than 0.1 px");
		check(sized && worst_kept_residual < 0.02f, "a kept subset's offset is reported as the fit's own small residual");
		const double check_apart = worstAgainst(r.projective, r.transform, ref_w, ref_h);
		printf("ground partly on the specimen: the perspective check departs %.4f px from the similarity\n", check_apart);
		check(r.projective_fitted && check_apart < 0.01,
			"the perspective check leaves out the ground that moved too, and agrees with the similarity");
	}

	// 12. Homography2D's own arithmetic, by hand. H maps (X, Y) to
	//     ((2X + 1) / w, (3Y - 2) / w) with w = 0.01 X + 1: at (10, 4), w = 1.1 and
	//     the point goes to (21 / 1.1, 10 / 1.1).
	{
		Homography2D H;
		const double h[9] = { 2, 0, 1, 0, 3, -2, 0.01, 0, 1 };
		for (int i = 0; i < 9; i++) H.h[i] = h[i];
		const Point2D m = H.map(Point2D(10, 4));
		check(fabs(m.x - 21 / 1.1) < 1e-4 && fabs(m.y - 10 / 1.1) < 1e-4, "homography map: (10, 4) goes to (19.09, 9.09)");
		const Point2D u = H.unmap(m);
		check(fabs(u.x - 10) < 1e-4 && fabs(u.y - 4) < 1e-4, "homography unmap undoes map");
		//d/dX of (2X + 1) / (0.01 X + 1) at 10 is (2 * 1.1 - 21 * 0.01) / 1.21
		double j[4];
		H.jacobian(Point2D(10, 4), j);
		check(fabs(j[0] - (2.2 - 0.21) / 1.21) < 1e-6 && fabs(j[1]) < 1e-9
			&& fabs(j[2] + 10 * 0.01 / 1.21) < 1e-6 && fabs(j[3] - 3 / 1.1) < 1e-6,
			"homography jacobian at (10, 4), worked by hand");

		//a POI seeded from H and measured as exactly H carries back to no movement at all
		POI2D poi(Point2D(30, 17));
		H.seed(poi);
		H.toReference(poi);
		check(fabs(poi.deformation.u) < 1e-4 && fabs(poi.deformation.v) < 1e-4 && fabs(poi.deformation.ux) < 1e-5
			&& fabs(poi.deformation.uy) < 1e-5 && fabs(poi.deformation.vx) < 1e-5 && fabs(poi.deformation.vy) < 1e-5,
			"a point that moved only with the camera carries back to no movement");

		//a similarity as a homography maps every point as the similarity does
		Similarity2D s;
		s.scale = 1.07; s.angle = 0.3; s.tx = 4; s.ty = -9;
		const Homography2D hs = Homography2D::from(s);
		const Point2D a = s.map(Point2D(57, 23)), b = hs.map(Point2D(57, 23));
		check(fabs(a.x - b.x) < 1e-4 && fabs(a.y - b.y) < 1e-4, "a similarity becomes the homography that maps as it does");
	}

	// 13. the camera TILTED between the sessions: the target's scale varies across
	//     the picture by about 5 %. Asked for perspective, the fit on the ground places
	//     every reference pixel within 0.01 px, and the specimen carried back through it
	//     reads its own deformation. The similarity on the same ground cannot: it is the
	//     reason the perspective fit exists.
	{
		TiltedCamera tilt;
		{
			const double c = cos(ANGLE), s = sin(ANGLE);
			const double base[9] = { SCALE * c, -SCALE * s, 0, SCALE * s, SCALE * c, 0, 2.0e-4, -1.5e-4, 1 };
			for (int i = 0; i < 9; i++) tilt.h[i] = base[i];
			const Point2D at = tilt.map(CX, CY);
			const double w = tilt.h[6] * CX + tilt.h[7] * CY + 1;
			tilt.h[2] = (140 - at.x) * w;
			tilt.h[5] = (120 - at.y) * w;
		}
		Image2D tilted_tar = render(tar_w, tar_h, nullptr, true, specimen, rig, &tilt);

		SessionAlignment2D flat(15, 10, 4);
		flat.setImages(ref, tilted_tar);
		flat.setGround(&ground);
		const AlignmentResult f = flat.compute();

		SessionAlignment2D persp(15, 10, 4);
		persp.setImages(ref, tilted_tar);
		persp.setGround(&ground);
		persp.setPerspective(true);
		const AlignmentResult p = persp.compute();

		const double worst_f = worstAgainstTilted(f.transform, tilt, ref_w, ref_h);
		const double worst_p = worstAgainstTilted(p.projective, tilt, ref_w, ref_h);
		printf("tilted, similarity: status %d, %d of %d used, rms %.4f px, worst mapping error %.4f px\n",
			(int)f.status, f.refined_used, f.refined_points, f.residual_rms, worst_f);
		printf("tilted, perspective: status %d, %d of %d used, rms %.4f px, worst mapping error %.5f px\n",
			(int)p.status, p.refined_used, p.refined_points, p.residual_rms, worst_p);
		check(p.status == AlignmentStatus::ok && p.perspective, "asked for perspective, the fit is a homography");
		check(!f.perspective, "not asked, it is not");
		check(worst_p < 0.01, "a tilted camera is recovered within 0.01 px anywhere on the image");
		check(worst_f > 0.3, "the similarity cannot describe a tilted camera");
		check(p.residual_rms < 0.01f && p.refined_used == p.refined_points,
			"the residual and counts reported are the perspective fit's, not the similarity's");
		//and unasked, the perspective check through the similarity's own subsets tells
		//the tilt apart: it departs from the similarity by pixels and sits near the truth,
		//which is what lets a caller say "the camera may have tilted" from one run
		const double check_apart = worstAgainst(f.projective, f.transform, ref_w, ref_h);
		const double check_truth = worstAgainstTilted(f.projective, tilt, ref_w, ref_h);
		printf("tilted, unasked: the perspective check departs %.3f px from the similarity, %.4f px from the truth\n",
			check_apart, check_truth);
		check(f.projective_fitted && check_apart > 0.3 && check_truth < 0.05,
			"unasked, the perspective check through the same subsets tells a tilted camera from a similarity");

		//the specimen through the perspective fit, measured with second-order subsets: under
		//a tilt the camera's warp across a subset curves, and first-order subsets measure
		//this specimen 0.023 px off on average for that reason alone
		ICGN2D2 icgn(16, 16, 0.001f, 20, 4);
		icgn.setImages(ref, tilted_tar);
		icgn.prepare();
		vector<POI2D> pois;
		for (int y = (int)CY - 36; y <= (int)CY + 36; y += 6)
			for (int x = (int)CX - 36; x <= (int)CX + 36; x += 6)
			{
				if ((x - CX) * (x - CX) + (y - CY) * (y - CY) > 40 * 40) continue;
				POI2D poi(Point2D((float)x, (float)y));
				p.projective.seed(poi);
				pois.push_back(poi);
			}
		icgn.compute(pois);
		double error = 0, g[4] = { 0, 0, 0, 0 };
		int solved = 0;
		for (auto& poi : pois)
		{
			if (poi.result.zncc < 0.9f) continue;
			p.projective.toReference(poi);
			const Point2D t = deform(poi.x, poi.y);
			error += hypot(poi.deformation.u - (t.x - poi.x), poi.deformation.v - (t.y - poi.y));
			g[0] += poi.deformation.ux; g[1] += poi.deformation.uy; g[2] += poi.deformation.vx; g[3] += poi.deformation.vy;
			solved++;
		}
		error /= max(solved, 1);
		for (double& k : g) k /= max(solved, 1);
		printf("tilted specimen: %d of %zu solved, mean displacement error %.5f px, gradient %+.5f %+.5f %+.5f %+.5f\n",
			solved, pois.size(), error, g[0], g[1], g[2], g[3]);
		check(solved == (int)pois.size() && error < 0.01, "through the perspective fit the specimen's displacement is within 0.01 px");
		check(fabs(g[0] - E11) < 3e-4 && fabs(g[1] - E12) < 3e-4 && fabs(g[2] - E21) < 3e-4 && fabs(g[3] - E22) < 3e-4,
			"and each gradient term within 3e-4 of E");
	}

	// 14. not asked for perspective, the homography through the same ground is still
	//     fitted, as a check on the similarity: on the untilted pair it agrees with the
	//     similarity everywhere, so a caller can say the camera did not tilt.
	{
		const double apart = worstAgainst(on_ground.projective, on_ground.transform, ref_w, ref_h);
		printf("untilted: the perspective check departs from the similarity by %.5f px at most\n", apart);
		check(on_ground.projective_fitted, "the perspective check is fitted alongside the similarity");
		check(apart < 0.01, "and on an untilted camera it agrees with the similarity within 0.01 px");
	}

	printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}
