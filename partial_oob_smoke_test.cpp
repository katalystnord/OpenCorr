/*
 Local build-verification test for the partially-out-of-bounds target subset
 guard (upstream PR #25, merged here 2026-08-17).

 Why this test exists, specifically: #25 was authored on a standalone branch
 off upstream and merged there, but never merged back into this fork, so
 surview-dev silently lacked the fix for two weeks. Nothing caught it, because
 no test exercised the guard. status_flag_smoke_test's own "out-of-bounds"
 case looks like it covers this but does not -- it places a POI too close to
 the image edge, which trips the *reference* subset bounds pre-check at the
 top of compute(). That check predates #25. The case #25 actually guards is
 different and strictly harder: a reference subset comfortably inside the
 image, whose *warped target* subset is pushed partially off the image by the
 incoming deformation guess.

 Interpolation2D::compute() reports an out-of-range coordinate by returning
 BicubicBspline::OUT_OF_BOUNDS (-1.f), not NaN. Unguarded, those fabricated
 -1 values are averaged into the subset's own mean and norm and blended into
 the ZNSSD cost as if they were real pixel intensities.

 What that actually costs, measured rather than assumed: with the guards
 temporarily removed, every partial-overlap case probed here (down to the
 mildest possible -- a single out-of-range column, 33 of 1089 samples)
 returns STATUS_MAX_ITERATIONS_REACHED, not a plausible ZNCC. A fabricated
 -1 sitting among real 0-255 intensities is a large enough outlier to keep
 the solver from converging at all. So the harm is not a convincing wrong
 number; it is a misdiagnosed one -- the POI is reported as a hard
 correlation failure ("did not converge") when the real cause is a
 deformation guess that pushed the subset off the image, which is a
 different problem with a different fix. The guard converts that into the
 accurate code. Worth keeping in mind if this is ever revisited: nothing
 guarantees the corrupted cost always fails to converge, and a milder
 overlap on a lower-contrast image could plausibly land on a real-looking
 ZNCC instead.

 Seven checks:
 1. The test's own geometry is really the partial case -- probing the
    interpolator directly over the exact sample grid ICGN2D1 will use, the
    warped target subset must contain BOTH valid samples AND out-of-bounds
    ones. Without this, a future change to the interpolator's border margin
    could quietly make the case fully in-bounds, leaving the rest passing
    while testing nothing at all.
 2-3. ICGN2D1: the partially-out-of-bounds guess is rejected as
    STATUS_INVALID_SUBSET_OR_GUESS, while the same POI with an in-bounds
    guess still correlates normally (so the guard is not simply firing on
    everything).
 4. The minimal case: a guess that puts just one column out of range must be
    rejected too, not only a grossly off-image one.
 5-7. The same rejection at the other three 2D guard sites: ICGN2D2, and the
    center_offset overloads of both (used by the stereo path).

 Not covered here: ICGN3D1 carries the same guard, but reaches it by a
 different mechanism -- a flag set during its sampling loop rather than a
 post-hoc Eigen .any() check, since a 3D subset is a vol_mat nested vector
 rather than an Eigen array. Exercising it needs a volume, not the 2D example
 images this test uses.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

//geometry shared by every check below, chosen against the 280x900 example
//images so that the reference subset sits well inside the image while the
//guess pushes the target subset off its left edge
const int SUBSET_RADIUS = 16;
const float POI_X = 20.f;   //ref subset spans x in [4, 36] -- comfortably in bounds
const float POI_Y = 400.f;  //ref subset spans y in [384, 416] -- comfortably in bounds
const float GUESS_U = -18.f;      //target subset spans x in [-14, 18] -- clearly off-image
const float GUESS_U_MINIMAL = -4.f; //target subset spans x in [0, 32] -- only x=0 unsamplable

int main()
{
	int failures = 0;

	Image2D ref_img("examples/2d_dic/oht_cfrp_0.bmp");
	Image2D tar_img("examples/2d_dic/oht_cfrp_4.bmp");

	//--- check 1: the geometry really is the *partial* case ---
	cout << "=== Geometry: warped target subset must be partially, not wholly, out of bounds ===" << endl;
	BicubicBspline probe(tar_img);
	probe.prepare();

	int oob_samples = 0, valid_samples = 0;
	for (int y_local = -SUBSET_RADIUS; y_local <= SUBSET_RADIUS; y_local++)
	{
		for (int x_local = -SUBSET_RADIUS; x_local <= SUBSET_RADIUS; x_local++)
		{
			//identity shape function apart from the u guess, matching what ICGN2D1
			//samples on its first iteration
			Point2D sample(POI_X + x_local + GUESS_U, POI_Y + y_local);
			if (probe.compute(sample) == BicubicBspline::OUT_OF_BOUNDS)
			{
				oob_samples++;
			}
			else
			{
				valid_samples++;
			}
		}
	}
	bool geometry_ok = oob_samples > 0 && valid_samples > 0;
	cout << "  " << valid_samples << " valid + " << oob_samples << " out-of-bounds sample(s) "
		<< (geometry_ok ? "PASS" : "FAIL -- case is degenerate, checks below prove nothing") << endl;
	if (!geometry_ok) failures++;

	//--- checks 2-3: ICGN2D1, plain compute() ---
	cout << endl << "=== ICGN2D1::compute(POI2D*) ===" << endl;
	ICGN2D1 icgn1(SUBSET_RADIUS, SUBSET_RADIUS, 0.001f, 20, 1);
	icgn1.setImages(ref_img, tar_img);
	icgn1.prepare();

	POI2D poi_partial(Point2D(POI_X, POI_Y));
	poi_partial.deformation.u = GUESS_U;
	icgn1.compute(&poi_partial);
	bool partial_ok = poi_partial.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  partially out-of-bounds guess -> zncc=" << poi_partial.result.zncc
		<< " (\"" << statusDescription(poi_partial.result.zncc) << "\") "
		<< (partial_ok ? "PASS" : "FAIL") << endl;
	if (!partial_ok) failures++;

	//control: identical POI and subset, in-bounds guess -- must still correlate
	POI2D poi_control(Point2D(POI_X, POI_Y));
	icgn1.compute(&poi_control);
	bool control_ok = !isFailureStatus(poi_control.result.zncc);
	cout << "  in-bounds guess, same POI    -> zncc=" << poi_control.result.zncc
		<< " (\"" << statusDescription(poi_control.result.zncc) << "\") "
		<< (control_ok ? "PASS" : "FAIL") << endl;
	if (!control_ok) failures++;

	//--- check 4: the mildest partial overlap the geometry allows ---
	//one column of the warped subset falls outside the interpolator's border
	//margin, 33 of 1089 samples -- the guard must fire here too, not just on a
	//grossly off-image guess
	POI2D poi_minimal(Point2D(POI_X, POI_Y));
	poi_minimal.deformation.u = GUESS_U_MINIMAL;
	icgn1.compute(&poi_minimal);
	bool minimal_ok = poi_minimal.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  single out-of-range column   -> zncc=" << poi_minimal.result.zncc
		<< " (\"" << statusDescription(poi_minimal.result.zncc) << "\") "
		<< (minimal_ok ? "PASS" : "FAIL") << endl;
	if (!minimal_ok) failures++;

	//--- check 5: ICGN2D2, plain compute() ---
	cout << endl << "=== ICGN2D2::compute(POI2D*) ===" << endl;
	ICGN2D2 icgn2(SUBSET_RADIUS, SUBSET_RADIUS, 0.001f, 20, 1);
	icgn2.setImages(ref_img, tar_img);
	icgn2.prepare();

	POI2D poi_partial2(Point2D(POI_X, POI_Y));
	poi_partial2.deformation.u = GUESS_U;
	icgn2.compute(&poi_partial2);
	bool partial2_ok = poi_partial2.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  partially out-of-bounds guess -> zncc=" << poi_partial2.result.zncc
		<< " (\"" << statusDescription(poi_partial2.result.zncc) << "\") "
		<< (partial2_ok ? "PASS" : "FAIL") << endl;
	if (!partial2_ok) failures++;

	//--- checks 6-7: the center_offset overloads (stereo path) ---
	//the offset cancels out of the sampled global coordinate under an identity
	//shape function (center moves by +offset, local coordinates by -offset), so
	//the same guess drives the subset off the image edge by the same amount
	cout << endl << "=== center_offset overloads ===" << endl;
	Point2D center_offset(5.f, 0.f);

	POI2D poi_offset1(Point2D(POI_X, POI_Y));
	poi_offset1.deformation.u = GUESS_U;
	icgn1.compute(&poi_offset1, center_offset);
	bool offset1_ok = poi_offset1.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  ICGN2D1 + center_offset -> zncc=" << poi_offset1.result.zncc
		<< " (\"" << statusDescription(poi_offset1.result.zncc) << "\") "
		<< (offset1_ok ? "PASS" : "FAIL") << endl;
	if (!offset1_ok) failures++;

	POI2D poi_offset2(Point2D(POI_X, POI_Y));
	poi_offset2.deformation.u = GUESS_U;
	icgn2.compute(&poi_offset2, center_offset);
	bool offset2_ok = poi_offset2.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  ICGN2D2 + center_offset -> zncc=" << poi_offset2.result.zncc
		<< " (\"" << statusDescription(poi_offset2.result.zncc) << "\") "
		<< (offset2_ok ? "PASS" : "FAIL") << endl;
	if (!offset2_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")
		<< " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
