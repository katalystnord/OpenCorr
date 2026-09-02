/*
 Local build-verification test for RegionFit2D, the Regional Fitting module
 upstream added on 2026.08.14 as a remedial procedure for POIs where ICGN
 fails from a poor initial guess: it fits an affine displacement field to the
 reliable POIs around such a point and writes the result back as that point's
 deformation.

 Written here for two reasons. First, the module arrived after both of this
 fork's full-branch audits (#17, #18), so nothing in this branch had ever
 exercised it. Second, and the reason it exists at all: RegionFit2D was
 unreachable through opencorr.h, the umbrella header every consumer of this
 library includes -- oc_region_fit.h was simply missing from it, so the class
 compiled into the library and could not be named by anyone using it the
 intended way. Case 1 is that regression, expressed as a test.

 1. Reachability. This file includes opencorr.h and nothing else, so it fails
    to COMPILE if the umbrella header stops exporting the module. Verified by
    removing the include again: "'RegionFit2D' is not a member of 'opencorr'".

 2. Exact recovery of an affine field. The reliable neighbours carry
    u = 2 + 0.1x - 0.05y and v = -1 + 0.02x + 0.08y exactly, so the least
    squares fit through them is that same plane with no residual, and the
    recovered u, ux, uy, v, vx, vy must match the coefficients to within
    float noise. This checks the local-coordinate convention too: the fit is
    solved about the POI's own position, so u is the value AT the POI and not
    at the origin.

 3. zncc is reset to 0 on success. The fitted deformation is an estimate
    borrowed from the neighbours, not a measurement of this point, and the
    module marks it as such. A caller that carried the old zncc forward would
    report an unmeasured point as well correlated.

 4. ⚑ A point it declines to fit is left EXACTLY as it was, silently. With
    fewer reliable POIs in total than neighbor_number_min, both the radius
    search and the KNN fallback come up short and compute() returns having
    touched nothing. "Declined" and "untouched" are therefore the same state,
    and a caller cannot tell them apart from the POI alone -- it must compare
    against what it passed in, or track the attempt itself. The same shape of
    trap as a rejected POI holding a displacement of zero.
*/

#include <cmath>
#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

namespace
{
	//the affine field the reliable neighbours are given, exactly
	const float u0 = 2.f, ux = 0.1f, uy = -0.05f;
	const float v0 = -1.f, vx = 0.02f, vy = 0.08f;

	float trueU(float x, float y) { return u0 + ux * x + uy * y; }
	float trueV(float x, float y) { return v0 + vx * x + vy * y; }

	bool close(float a, float b, float tol)
	{
		return fabs(a - b) <= tol;
	}

	int check(const char* what, bool ok)
	{
		cout << (ok ? "  OK   " : "  FAIL ") << what << endl;
		return ok ? 0 : 1;
	}
}

int main()
{
	int failures = 0;

	//--- a grid of reliable POIs carrying the affine field ---
	cout << "=== Regional fitting over reliable neighbours ===" << endl;
	vector<POI2D> reliable;
	for (int y = 80; y <= 120; y += 10)
	{
		for (int x = 80; x <= 120; x += 10)
		{
			POI2D poi(x, y);
			poi.deformation.u = trueU((float)x, (float)y);
			poi.deformation.v = trueV((float)x, (float)y);
			poi.result.zncc = 0.99f;
			reliable.push_back(poi);
		}
	}
	cout << "  " << reliable.size() << " reliable POIs" << endl;

	//radius wide enough to hold a good neighbourhood, minimum well below it
	RegionFit2D region_fit(25.f, 6, 1);
	region_fit.setNeighbor(reliable);
	region_fit.prepare();

	//--- 2 and 3: a point in the middle of them, as if ICGN had failed here ---
	POI2D failed(100, 100);
	failed.deformation.u = 0.f;
	failed.deformation.v = 0.f;
	failed.result.zncc = -2.f; //the sentinel a failed solve leaves behind
	region_fit.compute(&failed);

	//float tolerance: the fit is a QR solve over 25 points, not exact arithmetic
	const float tol = 1e-3f;
	failures += check("u recovered at the POI's own position",
		close(failed.deformation.u, trueU(100.f, 100.f), tol));
	failures += check("v recovered at the POI's own position",
		close(failed.deformation.v, trueV(100.f, 100.f), tol));
	failures += check("ux recovered", close(failed.deformation.ux, ux, tol));
	failures += check("uy recovered", close(failed.deformation.uy, uy, tol));
	failures += check("vx recovered", close(failed.deformation.vx, vx, tol));
	failures += check("vy recovered", close(failed.deformation.vy, vy, tol));
	failures += check("zncc reset to 0, marking the value as fitted not measured",
		failed.result.zncc == 0.f);

	cout << "  fitted u = " << failed.deformation.u
		<< ", expected " << trueU(100.f, 100.f) << endl;

	//--- 4: too few reliable POIs to fit at all ---
	cout << "=== A point it declines to fit ===" << endl;
	vector<POI2D> too_few;
	for (int i = 0; i < 3; i++)
	{
		POI2D poi(500 + i * 10, 500);
		poi.deformation.u = trueU(500.f + i * 10, 500.f);
		poi.deformation.v = trueV(500.f + i * 10, 500.f);
		too_few.push_back(poi);
	}

	//minimum of 6 against 3 available: neither the radius search nor the KNN
	//fallback can reach it
	RegionFit2D starved(25.f, 6, 1);
	starved.setNeighbor(too_few);
	starved.prepare();

	POI2D untouched(505, 500);
	untouched.deformation.u = 7.f;
	untouched.deformation.v = 8.f;
	untouched.result.zncc = -2.f;
	starved.compute(&untouched);

	failures += check("u left exactly as it was", untouched.deformation.u == 7.f);
	failures += check("v left exactly as it was", untouched.deformation.v == 8.f);
	failures += check("zncc left as it was, so a decline is silent",
		untouched.result.zncc == -2.f);
	cout << "  a declined point is indistinguishable from an untried one;" << endl;
	cout << "  the caller must track the attempt, not read it back." << endl;

	cout << endl;
	if (failures != 0)
	{
		cout << failures << " check(s) FAILED" << endl;
		return 1;
	}
	cout << "region_fit_smoke_test: all checks passed" << endl;
	return 0;
}
