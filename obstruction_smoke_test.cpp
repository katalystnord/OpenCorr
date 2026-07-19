/*
 Local build-verification test for Subset2D's obstruction/occlusion masking
 data model, issue #14.

 Phase 1 only (see oc_subset.h's scope note): this tests the mark/clear/query
 API itself, plus confirms marking pixels obstructed has NO effect yet on
 fill()/zeroMeanNorm() -- phase 2 (threading obstruction-awareness into the
 actual correlation math) is a separate, not-yet-scoped effort.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	//--- check 1: default state, mark/query, clear individual, clear all ---
	cout << "=== mark/clear/query API ===" << endl;
	Subset2D subset(Point2D(50.f, 50.f), 10, 10);

	bool default_clear = !subset.isObstructed(10, 10) && !subset.isObstructed(0, 0);
	cout << "  fresh subset has no obstructed pixels: " << (default_clear ? "PASS" : "FAIL") << endl;
	if (!default_clear) failures++;

	subset.markObstructed(42, 43);
	subset.markObstructed(44, 45);
	bool marked_ok = subset.isObstructed(42, 43) && subset.isObstructed(44, 45) && !subset.isObstructed(1, 1);
	cout << "  markObstructed() + isObstructed() roundtrip: " << (marked_ok ? "PASS" : "FAIL") << endl;
	if (!marked_ok) failures++;

	subset.clearObstructed(42, 43);
	bool clear_one_ok = !subset.isObstructed(42, 43) && subset.isObstructed(44, 45);
	cout << "  clearObstructed() removes only the given pixel: " << (clear_one_ok ? "PASS" : "FAIL") << endl;
	if (!clear_one_ok) failures++;

	subset.clearAllObstructed();
	bool clear_all_ok = !subset.isObstructed(44, 45);
	cout << "  clearAllObstructed() removes everything: " << (clear_all_ok ? "PASS" : "FAIL") << endl;
	if (!clear_all_ok) failures++;

	//--- check 2: phase 1 truly has no effect on fill()/zeroMeanNorm() ---
	cout << endl << "=== regression: marking pixels obstructed doesn't change fill()/zeroMeanNorm() yet ===" << endl;
	Image2D ref_img("examples/2d_dic/oht_cfrp_0.bmp");

	Subset2D plain(Point2D(100.f, 100.f), 16, 16);
	plain.fill(&ref_img);
	float plain_norm = plain.zeroMeanNorm();

	Subset2D obstructed(Point2D(100.f, 100.f), 16, 16);
	//mark most of the subset's own footprint obstructed -- if this silently changed
	//fill()/zeroMeanNorm() behavior, that would be phase 2 sneaking in unannounced
	for (int y = 84; y <= 116; y++)
		for (int x = 84; x <= 116; x++)
			obstructed.markObstructed(x, y);
	obstructed.fill(&ref_img);
	float obstructed_norm = obstructed.zeroMeanNorm();

	bool no_effect_yet = plain_norm == obstructed_norm;
	cout << "  plain norm=" << plain_norm << " vs obstructed-but-phase-1 norm=" << obstructed_norm
		<< " " << (no_effect_yet ? "PASS" : "FAIL") << endl;
	if (!no_effect_yet) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
