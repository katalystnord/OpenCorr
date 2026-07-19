/*
 Local build-verification test for Strain, which had zero existing test
 coverage in this fork despite being touched by the degenerate-fit guard
 added alongside issue #16 (crack/discontinuity diagnostic).

 1. Known-answer case: a globally linear (affine) displacement field is
    assigned to every POI directly (u = ux*x + uy*y, v = vx*x + vy*y for
    chosen constants). A local [1,dx,dy] least-squares fit of a globally
    linear field is exact everywhere (no approximation error), so
    Strain::compute() should recover the exact ux/uy/vx/vy at every POI --
    a much stronger check than a loose numerical tolerance.
 2. Degenerate-fit guard: a query POI whose neighbors are all placed exactly
    collinear (same y, varying x) makes the [1,dx,dy] design matrix's dy
    column identically zero -- rank-deficient, less than 3. Confirms the
    guard added to Strain::compute() actually fires (leaves strain at its
    untouched default) instead of solving into an arbitrary result.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;
	omp_set_num_threads(4);

	//--- known-answer case: globally linear displacement field ---
	cout << "=== Known-answer case: globally linear (affine) displacement field ===" << endl;
	float true_ux = 0.02f, true_uy = -0.01f, true_vx = 0.005f, true_vy = 0.03f;

	vector<POI2D> poi_queue;
	for (int y = 0; y <= 100; y += 10)
	{
		for (int x = 0; x <= 100; x += 10)
		{
			POI2D poi(Point2D((float)x, (float)y));
			poi.deformation.u = true_ux * x + true_uy * y;
			poi.deformation.v = true_vx * x + true_vy * y;
			poi.result.zncc = 1.f;
			poi_queue.push_back(poi);
		}
	}

	Strain strain(25.f, 4, 4);
	strain.setApproximation(1); //Cauchy strain: exx=ux, eyy=vy, exy=0.5*(uy+vx), exact for this affine field
	strain.prepare(poi_queue);
	strain.compute(poi_queue);

	float max_err_exx = 0.f, max_err_eyy = 0.f, max_err_exy = 0.f;
	int interior_checked = 0;
	for (auto& poi : poi_queue)
	{
		//skip the queue's own border POIs -- their neighbor count within the search radius
		//is legitimately smaller and less numerically exact near the edge of the sampled grid
		if (poi.x < 20.f || poi.x > 80.f || poi.y < 20.f || poi.y > 80.f) continue;

		float expected_exx = true_ux, expected_eyy = true_vy, expected_exy = 0.5f * (true_uy + true_vx);
		max_err_exx = max(max_err_exx, fabs(poi.strain.exx - expected_exx));
		max_err_eyy = max(max_err_eyy, fabs(poi.strain.eyy - expected_eyy));
		max_err_exy = max(max_err_exy, fabs(poi.strain.exy - expected_exy));
		interior_checked++;
	}
	cout << "  checked " << interior_checked << " interior POIs, max abs error: exx=" << max_err_exx
		<< " eyy=" << max_err_eyy << " exy=" << max_err_exy << endl;
	bool affine_ok = interior_checked > 0 && max_err_exx < 1e-4f && max_err_eyy < 1e-4f && max_err_exy < 1e-4f;
	cout << "  " << (affine_ok ? "PASS" : "FAIL") << ": exact recovery of a globally linear field's strain" << endl;
	if (!affine_ok) failures++;

	//--- degenerate-fit guard ---
	cout << endl << "=== Degenerate-fit guard: exactly collinear neighbors ===" << endl;
	vector<POI2D> degenerate_queue;
	POI2D query_poi(Point2D(50.f, 50.f));
	query_poi.result.zncc = 1.f;
	degenerate_queue.push_back(query_poi);
	//every neighbor at the same y as the query POI: the dy column of the design matrix is
	//identically zero, so the fit can never distinguish uy/vy -- rank 2, not 3
	for (int dx = -20; dx <= 20; dx += 5)
	{
		if (dx == 0) continue;
		POI2D neighbor(Point2D(50.f + dx, 50.f));
		neighbor.deformation.u = 0.01f * dx;
		neighbor.deformation.v = 0.01f * dx;
		neighbor.result.zncc = 1.f;
		degenerate_queue.push_back(neighbor);
	}

	Strain degenerate_strain(25.f, 4, 4);
	degenerate_strain.prepare(degenerate_queue);
	degenerate_strain.compute(degenerate_queue);

	bool guard_ok = degenerate_queue[0].strain.exx == 0.f
		&& degenerate_queue[0].strain.eyy == 0.f
		&& degenerate_queue[0].strain.exy == 0.f;
	cout << "  query POI strain after degenerate fit: exx=" << degenerate_queue[0].strain.exx
		<< " eyy=" << degenerate_queue[0].strain.eyy << " exy=" << degenerate_queue[0].strain.exy << endl;
	cout << "  " << (guard_ok ? "PASS" : "FAIL") << ": degenerate (rank-deficient) local fit left strain untouched, not solved into a value" << endl;
	if (!guard_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
