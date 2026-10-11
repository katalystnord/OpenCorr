/*
 Local build-verification test for masked subsets, issue #19 (phase 2 of #14):
 correlating only over the pixels a ValidPixelMask2D keeps.

 The fixture is rendered here, analytically, so its answer is exact by
 construction and nothing is resampled. A specimen speckle carries a 40 px hole
 that moves WITH the specimen; through the hole shows a second, static speckle,
 the background, which does not move at all. That is the case that does the
 damage: unmasked, a subset beside the hole sees part of a pattern that stayed
 put, and its answer is pulled toward no movement while it still correlates
 well.

 Lopsided on purpose. The translation differs in x and y and in sign, the
 affine field has unequal and cross terms, and the points beside the hole sit
 in all four directions at several distances, so a mask applied transposed, or
 to the wrong side of the subset, cannot come out right by symmetry.

 Checks:
 1. OrderedRegion2D: the last shape containing a pixel decides, an island
    inside a cut is back in, a pixel in no shape is out, and the box is the
    added shapes' box.
 2. ValidPixelMask2D: pull-in is a square (Chebyshev) neighbourhood, so the
    diagonal pixel at (pull_in, pull_in) from an invalid one is invalid and
    the one at pull_in + 1 along an axis is valid; outside the image is valid;
    counts, shares and weights agree with isValid() pixel by pixel.
 3. ICGN2D1 and ICGN2D2, both compute() overloads, under translation and under
    an affine field: beside the hole the masked error is at the far-field
    level, while the SAME points unmasked are wrong by far more -- the second
    half proves the fixture bites.
 4. A subset the mask does not touch is computed bit-identically to a run with
    no mask at all.
 5. The minimum share: refused below it with STATUS_TOO_FEW_VALID_PIXELS, kept
    exactly at it.
 6. Out-of-bounds target samples at MASKED pixels do not reject a point; the
    same point unmasked is rejected, which proves the case is the right one.
 7. ICLM2D1, ICLM2D2, NR2D1 and SimplexMatch2D refuse a point whose subset the
    mask reaches with STATUS_MASK_NOT_SUPPORTED, and solve one it does not
    reach exactly as without a mask.
 8. Both new statuses are described and count as failures.
 9. Uncertainty2D: sigma over the valid pixels only (larger beside the hole
    than the full subset would claim), unchanged where the mask does not reach.
 10. FFTCC2D: a masked starting guess finds the specimen's movement where the
    unmasked one locks onto the static background.
 Plus: the masked mean and norm on a lopsided 3 x 5 subset by hand, and a swap of
 what shows through the hole, which must change nothing beside it.

 Written red first (2026-10-11), against stubs: 35 checks failed. NEGATIVE
 CHECKS, 28 mutations, each caught by the case named unless said: the first
 shape deciding instead of the last, a cut ignored, the box grown by a cut
 (OrderedRegion2D); the pull-in one short, a cross instead of a square, none,
 outside the image invalid, the weight transposed (ValidPixelMask2D); the floor
 compared with <=, or in double, or not at all; the mask never applied; the
 gradient unweighted, the out-of-bounds test unmasked, the target norm unmasked,
 each in all four ICGN bodies; the masked mean over the whole subset, invalid
 pixels left at minus the mean; subsetIsMasked() always false; the setter
 dropping the floor; a status undescribed; Uncertainty2D never masked, its
 gradient unweighted along either axis, its probes rejected off the image, its
 target norm unmasked; FFTCC2D unweighted, its target unweighted, its mean over
 the whole subset.

 Five first SURVIVED, and each needed a fixture rather than a rule. With the
 target exactly as bright as the reference, a mean over the wrong pixels left
 the same offset in both and cancelled, so the target is now brighter and of
 more contrast. Nothing changed what lay behind the hole, so a target norm that
 counted it went unseen; the swap test does. Sigma takes the weaker axis, and an
 unmasked axis only gains energy, so it is never the weaker one on an even
 speckle; speckles stretched along each axis make each the weaker in turn. And a
 wrong FFTCC mean biases the guess toward no movement, which a high-contrast
 speckle outweighs; a faint one on bright paint does not.

 One survives, and cannot be caught: sending every subset down the masked path.
 With nothing masked, that path computes the same bits as the unmasked one, so
 the bit-identity check passes either way. The unmasked path is kept for speed
 and for a diff against upstream that leaves its arithmetic untouched.
*/

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

namespace
{
	int failures = 0;

	void check(const string& name, bool ok, const string& detail = "")
	{
		cout << "  " << name << ": " << (ok ? "PASS" : "FAIL");
		if (!detail.empty())
		{
			cout << "  (" << detail << ")";
		}
		cout << endl;
		if (!ok)
		{
			failures++;
		}
	}

	//An analytic speckle: one Gaussian dot per 6 px cell, jittered, so a pixel
	//only needs the dots of its neighbouring cells.
	struct Speckle
	{
		int cells_x, cells_y;
		float amplitude, offset;
		float radius_x, radius_y; //unequal for a speckle stretched along one axis
		vector<float> dot_x, dot_y;
		static constexpr float cell = 6.f;

		Speckle(int width, int height, unsigned seed, float amplitude, float offset,
			float radius_x = 2.2f, float radius_y = 2.2f)
			: amplitude(amplitude), offset(offset), radius_x(radius_x), radius_y(radius_y)
		{
			//a margin of cells on every side, so a point warped off the image still
			//lands on pattern
			cells_x = (int)(width / cell) + 8;
			cells_y = (int)(height / cell) + 8;
			mt19937 rng(seed);
			uniform_real_distribution<float> jitter(0.f, cell);
			for (int j = 0; j < cells_y; j++)
			{
				for (int i = 0; i < cells_x; i++)
				{
					dot_x.push_back((i - 4) * cell + jitter(rng));
					dot_y.push_back((j - 4) * cell + jitter(rng));
				}
			}
		}

		float at(float x, float y) const
		{
			const int ci = (int)floor(x / cell) + 4;
			const int cj = (int)floor(y / cell) + 4;
			float sum = 0.f;
			for (int j = cj - 3; j <= cj + 3; j++)
			{
				for (int i = ci - 3; i <= ci + 3; i++)
				{
					if (i < 0 || j < 0 || i >= cells_x || j >= cells_y)
					{
						continue;
					}
					const float dx = x - dot_x[j * cells_x + i];
					const float dy = y - dot_y[j * cells_x + i];
					sum += exp(-(dx * dx) / (radius_x * radius_x) - (dy * dy) / (radius_y * radius_y));
				}
			}
			return offset + amplitude * sum;
		}
	};

	const int kWidth = 240, kHeight = 200;
	const int kHoleX = 120, kHoleY = 100;
	const float kHoleRadius = 20.f;
	const int kRadius = 16;

	//x' = F X + t, with X the material (reference) position
	struct Motion
	{
		float f11, f12, f21, f22, tx, ty;

		Point2D displacement(float x, float y) const
		{
			return Point2D(f11 * x + f12 * y + tx - x, f21 * x + f22 * y + ty - y);
		}
	};

	bool inHole(float x, float y)
	{
		const float dx = x - kHoleX, dy = y - kHoleY;
		return dx * dx + dy * dy <= kHoleRadius * kHoleRadius;
	}

	//The reference, and a target under `motion`. A material point in the hole shows
	//the background at whatever image position it occupies; everything else shows the
	//specimen pattern at its material position.
	//The target is also brighter and of more contrast than the reference, as a real
	//target is after the light drifts. A correlation is meant to be blind to that, and a
	//mean taken over the wrong pixels is not: with the two images equally bright, a wrong
	//mean cancelled between reference and target and every check passed against it.
	const float kTargetGain = 1.3f, kTargetOffset = 25.f;

	void render(const Motion& motion, const Speckle& specimen, const Speckle& background,
		Image2D& ref, Image2D& tar)
	{
		const float det = motion.f11 * motion.f22 - motion.f12 * motion.f21;
		for (int y = 0; y < kHeight; y++)
		{
			for (int x = 0; x < kWidth; x++)
			{
				ref.eg_mat(y, x) = inHole((float)x, (float)y) ? background.at((float)x, (float)y)
					: specimen.at((float)x, (float)y);

				const float px = x - motion.tx, py = y - motion.ty;
				const float mx = (motion.f22 * px - motion.f12 * py) / det;
				const float my = (-motion.f21 * px + motion.f11 * py) / det;
				tar.eg_mat(y, x) = kTargetOffset + kTargetGain * (inHole(mx, my) ? background.at((float)x, (float)y)
					: specimen.at(mx, my));
			}
		}
	}

	//What a region tool hands over for "the whole image, with the hole cut out".
	unique_ptr<OrderedRegion2D> imageWithHole()
	{
		auto region = make_unique<OrderedRegion2D>();
		region->add(make_unique<Polygon2D>(vector<int>{ 0, kWidth - 1, kWidth - 1, 0 },
			vector<int>{ 0, 0, kHeight - 1, kHeight - 1 }));
		region->cut(make_unique<Circle2D>(kHoleX, kHoleY, kHoleRadius));
		return region;
	}

	//Points beside the hole, `gap` px from its edge, in four directions; each is
	//inside the image with its whole subset.
	vector<Point2D> besideTheHole()
	{
		vector<Point2D> points;
		const int gaps[] = { 3, 5, 8 };
		for (int gap : gaps)
		{
			const int d = (int)kHoleRadius + gap;
			points.push_back(Point2D((float)(kHoleX + d), (float)(kHoleY + 1)));
			points.push_back(Point2D((float)(kHoleX - d), (float)(kHoleY - 2)));
			points.push_back(Point2D((float)(kHoleX + 2), (float)(kHoleY + d)));
			points.push_back(Point2D((float)(kHoleX - 1), (float)(kHoleY - d)));
		}
		return points;
	}

	//Points the mask cannot reach: every subset pixel more than pull-in from the hole.
	vector<Point2D> farFromTheHole()
	{
		return { Point2D(30.f, 30.f), Point2D(205.f, 40.f), Point2D(40.f, 165.f), Point2D(200.f, 170.f) };
	}

	vector<POI2D> seeded(const vector<Point2D>& points, const Motion& motion)
	{
		vector<POI2D> queue;
		for (const Point2D& p : points)
		{
			POI2D poi(p);
			const Point2D u = motion.displacement(p.x, p.y);
			//a starting guess within a pixel, as an integer-pixel search would give
			poi.deformation.u = round(u.x);
			poi.deformation.v = round(u.y);
			queue.push_back(poi);
		}
		return queue;
	}

	//mean |measured - exact| displacement over the points that solved; the number
	//that solved goes in `solved`
	float meanError(const vector<POI2D>& queue, const Motion& motion, int& solved)
	{
		float sum = 0.f;
		solved = 0;
		for (const POI2D& poi : queue)
		{
			if (!(poi.result.zncc > 0.f))
			{
				continue;
			}
			const Point2D u = motion.displacement(poi.x, poi.y);
			sum += hypot(poi.deformation.u - u.x, poi.deformation.v - u.y);
			solved++;
		}
		return solved > 0 ? sum / solved : 1.0e9f;
	}

	string describe(float error, int solved, size_t of)
	{
		return to_string(error) + " px, " + to_string(solved) + "/" + to_string(of) + " solved";
	}

	//Runs `solver` over `points`, masked or not, through the plain or the
	//center_offset overload.
	template<typename Solver>
	vector<POI2D> run(Solver& solver, const vector<Point2D>& points, const Motion& motion,
		const ValidPixelMask2D* mask, bool with_offset)
	{
		solver.setValidPixels(mask);
		vector<POI2D> queue = seeded(points, motion);
		if (with_offset)
		{
			vector<Point2D> offsets(queue.size(), Point2D(0.f, 0.f));
			solver.compute(queue, offsets);
		}
		else
		{
			solver.compute(queue);
		}
		return queue;
	}

	template<typename Solver>
	void accuracyChecks(const string& name, Solver& solver, Image2D& ref, Image2D& tar,
		const Motion& motion, const ValidPixelMask2D& mask, bool with_offset, float masked_limit)
	{
		solver.setImages(ref, tar);
		solver.prepare();
		const vector<Point2D> near = besideTheHole();

		int solved_masked = 0, solved_plain = 0;
		const float masked = meanError(run(solver, near, motion, &mask, with_offset), motion, solved_masked);
		const float plain = meanError(run(solver, near, motion, nullptr, with_offset), motion, solved_plain);

		check(name + ": masked, beside the hole, at the far-field level", masked < masked_limit
			&& solved_masked == (int)near.size(), describe(masked, solved_masked, near.size()));
		check(name + ": unmasked, the same points are pulled off (the fixture bites)",
			plain > 10.f * masked_limit, describe(plain, solved_plain, near.size()));

		//4. untouched subsets are bit-identical
		const vector<Point2D> far = farFromTheHole();
		const vector<POI2D> far_masked = run(solver, far, motion, &mask, with_offset);
		const vector<POI2D> far_plain = run(solver, far, motion, nullptr, with_offset);
		bool identical = true;
		for (size_t i = 0; i < far.size(); i++)
		{
			identical = identical && far_masked[i].deformation.u == far_plain[i].deformation.u
				&& far_masked[i].deformation.v == far_plain[i].deformation.v
				&& far_masked[i].deformation.ux == far_plain[i].deformation.ux
				&& far_masked[i].result.zncc == far_plain[i].result.zncc
				&& far_masked[i].result.zncc > 0.f;
		}
		check(name + ": a subset the mask does not touch is bit-identical", identical);
	}
}

int main()
{
	const Speckle specimen(kWidth, kHeight, 7u, 180.f, 30.f);
	const Speckle background(kWidth, kHeight, 1234u, 180.f, 30.f);

	//--- 1. OrderedRegion2D ---
	cout << "=== OrderedRegion2D: the last shape containing a pixel decides ===" << endl;
	{
		OrderedRegion2D region;
		region.add(make_unique<Polygon2D>(vector<int>{ 10, 90, 90, 10 }, vector<int>{ 20, 20, 60, 60 }));
		region.cut(make_unique<Circle2D>(50, 40, 15.f));
		region.add(make_unique<Circle2D>(52, 41, 4.f)); //an island back inside the cut
		region.cut(make_unique<Circle2D>(88, 58, 6.f)); //a cut reaching past the added shapes
		check("inside the first shape only: in", region.contains(15, 25));
		check("inside the cut: out", !region.contains(40, 40));
		check("inside the island added after the cut: in", region.contains(53, 42));
		check("in no shape: out", !region.contains(5, 5) && !region.contains(95, 40));
		check("four shapes listed", region.numShapes() == 4);
		check("the cut that reaches past them removes", !region.contains(87, 57));
		check("the box is the added shapes' box, not grown by a cut", region.getMinX() == 10 && region.getMaxX() == 90
			&& region.getMinY() == 20 && region.getMaxY() == 60);

		OrderedRegion2D cut_first;
		cut_first.cut(make_unique<Circle2D>(50, 40, 15.f));
		cut_first.add(make_unique<Polygon2D>(vector<int>{ 10, 90, 90, 10 }, vector<int>{ 20, 20, 60, 60 }));
		check("order matters: an addition after a cut covers it", cut_first.contains(40, 40));

		OrderedRegion2D empty;
		check("empty: contains nothing, and its box is empty",
			!empty.contains(0, 0) && empty.getMinX() > empty.getMaxX());
	}

	//--- 2. ValidPixelMask2D ---
	cout << endl << "=== ValidPixelMask2D: pull-in, edges, counts ===" << endl;
	{
		//one invalid pixel at (30, 20), in a lopsided image
		OrderedRegion2D region;
		region.add(make_unique<Polygon2D>(vector<int>{ 0, 69, 69, 0 }, vector<int>{ 0, 0, 39, 39 }));
		region.cut(make_unique<Circle2D>(30, 20, 0.5f)); //the one pixel (30, 20)
		const ValidPixelMask2D mask(region, 70, 40, 2);

		check("the cut pixel is invalid", !mask.isValid(30, 20));
		check("pulled in 2 px along x, both sides", !mask.isValid(32, 20) && !mask.isValid(28, 20)
			&& mask.isValid(33, 20) && mask.isValid(27, 20));
		check("pulled in 2 px along y, both sides", !mask.isValid(30, 22) && !mask.isValid(30, 18)
			&& mask.isValid(30, 23) && mask.isValid(30, 17));
		check("a square neighbourhood: the diagonal corner (2, 2) away is invalid",
			!mask.isValid(32, 22) && !mask.isValid(28, 18) && mask.isValid(33, 23));
		check("25 pixels invalid in all", mask.invalidCount() == 25, to_string(mask.invalidCount()));
		check("outside the image is valid", mask.isValid(-1, 5) && mask.isValid(70, 5) && mask.isValid(5, 40));

		//a 7 x 5 subset (radius 3, 2) centred at (34, 21): columns 31..37 and rows 19..23
		//take invalid columns 31, 32 over rows 19..22, which is 8 pixels
		check("valid count over a lopsided subset", mask.validCount(34, 21, 3, 2) == 35 - 8,
			to_string(mask.validCount(34, 21, 3, 2)));
		check("share", fabs(mask.validShare(34, 21, 3, 2) - 27.f / 35.f) < 1e-6f);
		Eigen::MatrixXf weight;
		const int n = mask.fillWeight(34, 21, 3, 2, weight);
		bool agrees = weight.rows() == 5 && weight.cols() == 7 && n == 27;
		for (int r = 0; agrees && r < 5; r++)
		{
			for (int c = 0; c < 7; c++)
			{
				agrees = agrees && (weight(r, c) == 1.f) == mask.isValid(34 + c - 3, 21 + r - 2);
			}
		}
		check("the weight is isValid(), rows by y", agrees);

		const ValidPixelMask2D no_pull(region, 70, 40, 0);
		check("no pull-in: only the cut pixel", no_pull.invalidCount() == 1 && no_pull.isValid(31, 20));
	}

	//--- 2b. the masked mean and norm ---
	cout << endl << "=== Subset2D: mean and norm over the valid pixels ===" << endl;
	{
		//a 3 x 5 subset, lopsided, with three pixels masked; the valid ones are
		//10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, mean 65
		Subset2D subset(Point2D(2.f, 1.f), 2, 1);
		subset.eg_mat << 10, 20, 999, 30, 40,
			50, 60, 70, -999, 80,
			90, 100, 110, 120, 555;
		Eigen::MatrixXf weight(3, 5);
		weight << 1, 1, 0, 1, 1,
			1, 1, 1, 0, 1,
			1, 1, 1, 1, 0;
		const float norm = subset.maskedZeroMeanNorm(weight, 12);
		float expected = 0.f;
		for (int v = 10; v <= 120; v += 10) expected += (v - 65.f) * (v - 65.f);
		check("the norm of the valid pixels about their own mean", fabs(norm - sqrt(expected)) < 1e-3f,
			to_string(norm) + ", " + to_string(sqrt(expected)) + " expected");
		check("valid pixels hold their difference from that mean", subset.eg_mat(0, 0) == -55.f
			&& subset.eg_mat(2, 3) == 55.f);
		check("masked pixels hold zero", subset.eg_mat(0, 2) == 0.f && subset.eg_mat(1, 3) == 0.f
			&& subset.eg_mat(2, 4) == 0.f);
	}

	//--- 3, 4. accuracy ---
	Image2D ref(kWidth, kHeight), tar_shift(kWidth, kHeight), tar_affine(kWidth, kHeight);
	const Motion shift{ 1.f, 0.f, 0.f, 1.f, 2.4f, -1.3f };
	const Motion affine{ 1.02f, 0.006f, -0.004f, 0.985f, 1.6f, 0.7f };
	render(shift, specimen, background, ref, tar_shift);
	{
		Image2D ref_again(kWidth, kHeight);
		render(affine, specimen, background, ref_again, tar_affine);
	}
	const auto region = imageWithHole();
	const ValidPixelMask2D mask(*region, kWidth, kHeight, 2);

	cout << endl << "=== ICGN2D1 ===" << endl;
	{
		ICGN2D1 icgn(kRadius, kRadius, 0.001f, 20, 1);
		accuracyChecks("translation", icgn, ref, tar_shift, shift, mask, false, 0.02f);
		accuracyChecks("translation, center_offset overload", icgn, ref, tar_shift, shift, mask, true, 0.02f);
		accuracyChecks("affine", icgn, ref, tar_affine, affine, mask, false, 0.03f);
	}
	cout << endl << "=== ICGN2D2 ===" << endl;
	{
		ICGN2D2 icgn(kRadius, kRadius, 0.001f, 20, 1);
		accuracyChecks("translation", icgn, ref, tar_shift, shift, mask, false, 0.02f);
		accuracyChecks("translation, center_offset overload", icgn, ref, tar_shift, shift, mask, true, 0.02f);
		accuracyChecks("affine", icgn, ref, tar_affine, affine, mask, false, 0.03f);
	}

	//--- 4b. what shows through the hole does not count ---
	cout << endl << "=== what shows through the hole takes no part ===" << endl;
	{
		//the same specimen and motion, a different background behind the hole, of another
		//brightness and contrast; masked, every figure beside the hole must stay put. To the
		//bit, measured: nothing a masked subset reads reaches into the hole once pulled in.
		//The one exception is sigma's NOISE estimate, which is taken once over the whole
		//reference image and so does see the hole (0.5 per cent here, 3 on a stretched
		//speckle); it is divided out, leaving the subset's own gradient energy.
		const Speckle other(kWidth, kHeight, 4321u, 120.f, 60.f);
		Image2D ref_other(kWidth, kHeight), tar_other(kWidth, kHeight);
		render(shift, specimen, other, ref_other, tar_other);

		const vector<Point2D> near = besideTheHole();
		auto measure = [&](Image2D& r, Image2D& t) {
			ICGN2D1 icgn(kRadius, kRadius, 0.001f, 20, 1);
			icgn.setImages(r, t);
			icgn.prepare();
			vector<POI2D> queue = run(icgn, near, shift, &mask, false);
			Uncertainty2D uncertainty(kRadius, kRadius, 1);
			uncertainty.setImages(r, t);
			uncertainty.prepare();
			uncertainty.setValidPixels(&mask);
			uncertainty.compute(queue);
			return queue;
		};
		const vector<POI2D> a = measure(ref, tar_shift);
		const vector<POI2D> b = measure(ref_other, tar_other);
		//the noise is estimated over the whole image, hole included, so it is divided out
		const float noise_a = Uncertainty2D::noiseStdDev(ref), noise_b = Uncertainty2D::noiseStdDev(ref_other);
		float du = 0.f, dz = 0.f, ds = 0.f, db = 0.f;
		for (size_t i = 0; i < near.size(); i++)
		{
			du = max(du, hypot(a[i].deformation.u - b[i].deformation.u, a[i].deformation.v - b[i].deformation.v));
			dz = max(dz, fabs(a[i].result.zncc - b[i].result.zncc));
			ds = max(ds, fabs((a[i].result.sigma / noise_a) / (b[i].result.sigma / noise_b) - 1.f));
			db = max(db, fabs(a[i].result.beta / b[i].result.beta - 1.f));
		}
		check("displacement, to the bit", du == 0.f, to_string(du) + " px apart");
		check("correlation, to the bit", dz == 0.f, to_string(dz) + " apart");
		check("match conditioning, to the bit", db == 0.f, to_string(db) + " apart");
		check("noise floor, once the whole-image noise estimate is divided out", ds < 1.0e-4f,
			to_string(100.f * ds) + " per cent apart");
	}

	//sigma takes the WEAKER of the two axes' gradient energy, and a masked pixel left in
	//can only add energy, so a mask applied to one axis alone is invisible on an even
	//speckle: the axis it inflates is never the weaker one. A speckle stretched along each
	//axis in turn makes each axis the weaker one.
	cout << endl << "=== the noise floor, on a speckle weaker along each axis ===" << endl;
	for (int axis = 0; axis < 2; axis++)
	{
		const Speckle stretched(kWidth, kHeight, 77u, 180.f, 30.f, axis == 0 ? 4.f : 1.6f, axis == 0 ? 1.6f : 4.f);
		const Speckle behind_a(kWidth, kHeight, 5u, 180.f, 30.f);
		const Speckle behind_b(kWidth, kHeight, 6u, 180.f, 30.f);
		Image2D ref_a(kWidth, kHeight), tar_a(kWidth, kHeight), ref_b(kWidth, kHeight), tar_b(kWidth, kHeight);
		render(shift, stretched, behind_a, ref_a, tar_a);
		render(shift, stretched, behind_b, ref_b, tar_b);
		auto sigma = [&](Image2D& r, Image2D& t) {
			Uncertainty2D uncertainty(kRadius, kRadius, 1);
			uncertainty.setImages(r, t);
			uncertainty.prepare();
			uncertainty.setValidPixels(&mask);
			vector<POI2D> queue = seeded(besideTheHole(), shift);
			for (POI2D& poi : queue) poi.result.zncc = 1.f; //as if solved at the seed
			uncertainty.compute(queue);
			return queue;
		};
		const vector<POI2D> a = sigma(ref_a, tar_a), b = sigma(ref_b, tar_b);
		//sigma = sqrt(2 noise^2 / energy), and the noise is estimated over the whole image,
		//hole included, so it is divided out: what is left is the gradient energy alone
		const float noise_a = Uncertainty2D::noiseStdDev(ref_a), noise_b = Uncertainty2D::noiseStdDev(ref_b);
		float ds = 0.f;
		for (size_t i = 0; i < a.size(); i++)
		{
			ds = max(ds, fabs((a[i].result.sigma / noise_a) / (b[i].result.sigma / noise_b) - 1.f));
		}
		check(string("weaker along ") + (axis == 0 ? "y" : "x") + ": the gradient energy does not see the hole",
			ds < 1.0e-4f, to_string(100.f * ds) + " per cent apart");
	}

	//--- 5. the minimum share ---
	cout << endl << "=== the minimum valid share ===" << endl;
	{
		ICGN2D1 icgn(kRadius, kRadius, 0.001f, 20, 1);
		icgn.setImages(ref, tar_shift);
		icgn.prepare();
		const Point2D p((float)(kHoleX + (int)kHoleRadius + 3), (float)(kHoleY + 1));
		const float share = mask.validShare((int)p.x, (int)p.y, kRadius, kRadius);
		check("the point is partly masked", share > 0.3f && share < 0.95f, to_string(share));

		vector<POI2D> refused = seeded({ p }, shift);
		icgn.setValidPixels(&mask, share + 0.01f);
		icgn.compute(refused);
		check("below the floor: STATUS_TOO_FEW_VALID_PIXELS",
			refused[0].result.zncc == (float)STATUS_TOO_FEW_VALID_PIXELS, to_string(refused[0].result.zncc));

		vector<POI2D> kept = seeded({ p }, shift);
		icgn.setValidPixels(&mask, share);
		icgn.compute(kept);
		check("exactly at the floor: kept", kept[0].result.zncc > 0.f, to_string(kept[0].result.zncc));

		ICGN2D1 fresh(kRadius, kRadius, 0.001f, 20, 1);
		check("the default floor is 25 %", fresh.min_valid_share == 0.25f
			&& DIC::kDefaultMinValidShare == 0.25f);
		fresh.setValidPixels(&mask);
		check("and setting a mask without a floor keeps it", fresh.min_valid_share == 0.25f);
		icgn.setValidPixels(nullptr);
		check("clearing the mask", icgn.valid_pixels == nullptr);
	}

	//--- 6. masked pixels off the target do not reject a point ---
	cout << endl << "=== out-of-bounds samples at masked pixels ===" << endl;
	{
		//cut everything right of x = 232, so the subset at x = 223 (columns 207..239)
		//keeps columns 207..229 after pull-in; the guess of +2 then carries the masked
		//columns off the right of the target, and the valid ones stay on it
		OrderedRegion2D right_cut;
		right_cut.add(make_unique<Polygon2D>(vector<int>{ 0, kWidth - 1, kWidth - 1, 0 },
			vector<int>{ 0, 0, kHeight - 1, kHeight - 1 }));
		right_cut.cut(make_unique<Polygon2D>(vector<int>{ 232, kWidth + 10, kWidth + 10, 232 },
			vector<int>{ -10, -10, kHeight + 10, kHeight + 10 }));
		const ValidPixelMask2D edge_mask(right_cut, kWidth, kHeight, 2);
		const Point2D p((float)(kWidth - 1 - kRadius), 60.f);

		ICGN2D1 icgn1(kRadius, kRadius, 0.001f, 20, 1);
		ICGN2D2 icgn2(kRadius, kRadius, 0.001f, 20, 1);
		icgn1.setImages(ref, tar_shift);
		icgn1.prepare();
		icgn2.setImages(ref, tar_shift);
		icgn2.prepare();

		int solved = 0;
		vector<POI2D> plain = run(icgn1, { p }, shift, nullptr, false);
		check("unmasked, the point is rejected for leaving the image (the case is the right one)",
			plain[0].result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS, to_string(plain[0].result.zncc));
		const float e1 = meanError(run(icgn1, { p }, shift, &edge_mask, false), shift, solved);
		check("ICGN2D1 masked: measured, accurately", solved == 1 && e1 < 0.02f, describe(e1, solved, 1));
		const float e1o = meanError(run(icgn1, { p }, shift, &edge_mask, true), shift, solved);
		check("ICGN2D1 center_offset masked: measured", solved == 1 && e1o < 0.02f, describe(e1o, solved, 1));
		const float e2 = meanError(run(icgn2, { p }, shift, &edge_mask, false), shift, solved);
		check("ICGN2D2 masked: measured, accurately", solved == 1 && e2 < 0.02f, describe(e2, solved, 1));
		const float e2o = meanError(run(icgn2, { p }, shift, &edge_mask, true), shift, solved);
		check("ICGN2D2 center_offset masked: measured", solved == 1 && e2o < 0.02f, describe(e2o, solved, 1));

		//and its reliability figures, whose probes carry the same masked pixels off the image
		vector<POI2D> measured = run(icgn1, { p }, shift, &edge_mask, false);
		Uncertainty2D uncertainty(kRadius, kRadius, 1);
		uncertainty.setImages(ref, tar_shift);
		uncertainty.prepare();
		uncertainty.setValidPixels(&edge_mask);
		uncertainty.compute(measured);
		check("Uncertainty2D masked: sigma and beta are measured there",
			measured[0].result.sigma > 0.f && measured[0].result.beta > 0.f,
			to_string(measured[0].result.sigma) + ", " + to_string(measured[0].result.beta));
	}

	//--- 7. solvers that do not honour a mask yet ---
	cout << endl << "=== solvers without mask support refuse, and only where it matters ===" << endl;
	{
		const vector<Point2D> near = { besideTheHole()[0] };
		const vector<Point2D> far = { farFromTheHole()[0] };

		auto refuses = [&](const string& name, auto& solver) {
			solver.setImages(ref, tar_shift);
			solver.prepare();
			solver.setValidPixels(&mask);
			vector<POI2D> touched = seeded(near, shift);
			solver.compute(touched);
			check(name + ": a masked subset is refused, saying why",
				touched[0].result.zncc == (float)STATUS_MASK_NOT_SUPPORTED, to_string(touched[0].result.zncc));

			vector<POI2D> untouched = seeded(far, shift);
			solver.compute(untouched);
			solver.setValidPixels(nullptr);
			vector<POI2D> without = seeded(far, shift);
			solver.compute(without);
			check(name + ": a subset the mask does not reach solves as without one",
				untouched[0].result.zncc > 0.f && untouched[0].deformation.u == without[0].deformation.u
				&& untouched[0].deformation.v == without[0].deformation.v,
				to_string(untouched[0].result.zncc) + " / " + to_string(without[0].result.zncc) + ", u "
				+ to_string(untouched[0].deformation.u) + " / " + to_string(without[0].deformation.u));
		};
		ICLM2D1 iclm1(kRadius, kRadius, 0.001f, 20, 1);
		refuses("ICLM2D1", iclm1);
		ICLM2D2 iclm2(kRadius, kRadius, 0.001f, 20, 1);
		refuses("ICLM2D2", iclm2);
		NR2D1 nr(kRadius, kRadius, 0.001f, 20, 1);
		refuses("NR2D1", nr);
		SimplexMatch2D simplex(kRadius, kRadius, 1);
		simplex.setIteration(2000, 1.0e-6f); //its default budget stops short on this fixture, mask or none
		refuses("SimplexMatch2D", simplex);
	}

	//--- 8. statuses ---
	cout << endl << "=== the new statuses ===" << endl;
	{
		const string unknown = statusDescription(-99.f);
		const string few = statusDescription((float)STATUS_TOO_FEW_VALID_PIXELS);
		const string unsupported = statusDescription((float)STATUS_MASK_NOT_SUPPORTED);
		check("too few valid pixels is described", few != unknown && few.find("valid") != string::npos, few);
		check("mask not supported is described", unsupported != unknown && unsupported != few, unsupported);
		check("both count as failures", isFailureStatus((float)STATUS_TOO_FEW_VALID_PIXELS)
			&& isFailureStatus((float)STATUS_MASK_NOT_SUPPORTED));
	}

	//--- 9. Uncertainty2D ---
	cout << endl << "=== Uncertainty2D: the noise floor of the pixels that took part ===" << endl;
	{
		ICGN2D1 icgn(kRadius, kRadius, 0.001f, 20, 1);
		icgn.setImages(ref, tar_shift);
		icgn.prepare();
		const Point2D near_point = besideTheHole()[0];
		const Point2D far_point = farFromTheHole()[0];
		vector<POI2D> solved = run(icgn, { near_point, far_point }, shift, &mask, false);

		Uncertainty2D uncertainty(kRadius, kRadius, 1);
		uncertainty.setImages(ref, tar_shift);
		uncertainty.prepare();

		vector<POI2D> masked = solved, plain = solved;
		uncertainty.setValidPixels(&mask);
		uncertainty.compute(masked);
		uncertainty.setValidPixels(nullptr);
		uncertainty.compute(plain);

		//sigma goes as 1/sqrt(gradient energy), and the energy of a share s of an even
		//speckle is about s times the whole: so masked/plain ~ 1/sqrt(s)
		const float share = mask.validShare((int)near_point.x, (int)near_point.y, kRadius, kRadius);
		const float ratio = masked[0].result.sigma / plain[0].result.sigma;
		const float expected = 1.f / sqrt(share);
		//the tolerance is the speckle's unevenness; a first version at 25 % passed against
		//an unmasked sigma, since 1/sqrt(0.58) is only 1.31, so the ratio must also clear 1.1
		check("beside the hole, sigma counts only the valid pixels",
			ratio > 1.1f && fabs(ratio / expected - 1.f) < 0.12f,
			"ratio " + to_string(ratio) + ", about " + to_string(expected) + " expected");
		check("beside the hole, beta is measured", masked[0].result.beta > 0.f, to_string(masked[0].result.beta));
		check("far from it, sigma and beta are unchanged", masked[1].result.sigma == plain[1].result.sigma
			&& masked[1].result.beta == plain[1].result.beta && masked[1].result.sigma > 0.f);
	}

	//--- 10. FFTCC2D ---
	cout << endl << "=== FFTCC2D: a masked starting guess ===" << endl;
	{
		//a large movement and a background of twice the contrast, so that unmasked the
		//static background outweighs the moving specimen in the subset
		const Motion big{ 1.f, 0.f, 0.f, 1.f, 7.f, -5.f };
		//and a faint specimen speckle on bright paint, which is where a mean taken over the
		//wrong pixels does its damage: it leaves the valid pixels offset by a constant,
		//and two offsets correlate best at no movement at all
		const Speckle loud(kWidth, kHeight, 99u, 240.f, 10.f);
		const Speckle faint(kWidth, kHeight, 7u, 50.f, 150.f);
		Image2D ref_big(kWidth, kHeight), tar_big(kWidth, kHeight);
		render(big, faint, loud, ref_big, tar_big);

		FFTCC2D fftcc(kRadius, kRadius, 1);
		fftcc.setImages(ref_big, tar_big);
		fftcc.prepare();
		const Point2D p((float)(kHoleX + (int)kHoleRadius + 5), (float)(kHoleY + 1));

		vector<POI2D> plain = { POI2D(p) };
		fftcc.compute(plain);
		vector<POI2D> masked = { POI2D(p) };
		fftcc.setValidPixels(&mask);
		fftcc.compute(masked);

		const string got_plain = to_string(plain[0].deformation.u) + ", " + to_string(plain[0].deformation.v);
		const string got_masked = to_string(masked[0].deformation.u) + ", " + to_string(masked[0].deformation.v);
		check("unmasked, the guess is pulled to the background (the fixture bites)",
			!(plain[0].deformation.u == 7.f && plain[0].deformation.v == -5.f), got_plain);
		check("masked, the guess is the specimen's movement",
			masked[0].deformation.u == 7.f && masked[0].deformation.v == -5.f, got_masked);
	}

	cout << endl << (failures == 0 ? "ALL PASS" : to_string(failures) + " FAILED") << endl;
	return failures == 0 ? 0 : 1;
}
