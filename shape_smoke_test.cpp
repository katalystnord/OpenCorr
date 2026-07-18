/*
 Local build-verification test for Shape2D/Polygon2D/Circle2D, issue #9.

 1. Known-geometry spot checks: a square (convex) and an L-shape (concave --
    the harder case for a point-in-polygon test) against hand-verified
    inside/outside/boundary points.
 2. Area cross-check: getOwnedPixels() pixel counts compared against the
    analytic area (shoelace formula for the polygons, pi*r^2 for the
    circle) -- a much stronger check than spot-testing a few points, since
    a systematically-biased contains() would still pass individual spot
    checks but fail an area comparison.
 3. Practical use: build a masked POI queue directly from getOwnedPixels(),
    confirming it's already usable for ROI-restricted POI construction even
    without kernel (per-subset) mask integration.
*/

#include <cmath>
#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	//--- square (convex), 10x10, corners at (10,10) and (20,20) ---
	cout << "=== Square (convex) ===" << endl;
	vector<int> sq_x = { 10, 20, 20, 10 };
	vector<int> sq_y = { 10, 10, 20, 20 };
	Polygon2D square(sq_x, sq_y);

	bool square_ok = square.contains(15, 15)   //center: inside
		&& square.contains(11, 11)             //just inside corner
		&& !square.contains(5, 15)             //left of the square: outside
		&& !square.contains(25, 15)            //right of the square: outside
		&& !square.contains(15, 5)             //above: outside
		&& !square.contains(15, 25);           //below: outside
	cout << "  " << (square_ok ? "PASS" : "FAIL") << ": spot checks (center in, well outside each direction out)" << endl;
	if (!square_ok) failures++;

	//--- L-shape (concave) ---
	// (0,0)-(10,0)-(10,5)-(5,5)-(5,10)-(0,10), an L covering the top-left 10x10 minus the
	// bottom-right 5x5 quadrant
	cout << endl << "=== L-shape (concave) ===" << endl;
	vector<int> l_x = { 0, 10, 10, 5, 5, 0 };
	vector<int> l_y = { 0, 0, 5, 5, 10, 10 };
	Polygon2D lshape(l_x, l_y);

	bool concave_ok = lshape.contains(2, 2)     //top-left arm: inside
		&& lshape.contains(8, 2)                //top-right arm: inside
		&& lshape.contains(2, 8)                //bottom-left arm: inside
		&& !lshape.contains(8, 8);              //the notched-out bottom-right quadrant: outside
	cout << "  " << (concave_ok ? "PASS" : "FAIL") << ": correctly handles the concave notch (not just a convex hull)" << endl;
	if (!concave_ok) failures++;

	//--- area cross-check ---
	cout << endl << "=== Area cross-check (pixel count vs analytic area) ===" << endl;
	auto shoelaceArea = [](const vector<int>& xs, const vector<int>& ys) -> float
	{
		float area = 0.f;
		int n = (int)xs.size();
		for (int i = 0; i < n; i++)
		{
			int j = (i + 1) % n;
			area += (float)xs[i] * ys[j] - (float)xs[j] * ys[i];
		}
		return fabs(area) / 2.f;
	};

	//Polygon2D now classifies every boundary pixel as inside (a deliberate fix for a prior
	//inconsistency where different edges of the same simple polygon resolved oppositely --
	//see oc_shape.cpp), so the discrete pixel count is the CLOSED-boundary count, not the
	//continuous shoelace area: an NxN axis-aligned square with all four edges included is
	//(N+1)^2 pixels, not N^2 -- the classic fencepost distinction between a continuous area
	//and an inclusive discrete grid. Asserting the exact hand-derived count (rather than a
	//loose tolerance around the continuous area) both confirms the fix and pins the
	//convention down precisely.
	float square_analytic = shoelaceArea(sq_x, sq_y);
	int square_pixels = (int)square.getOwnedPixels().size();
	int square_expected = 11 * 11; //vertices at (10,10)-(20,20): 11 integer values per axis, inclusive
	cout << "  square: analytic area=" << square_analytic << ", pixel count=" << square_pixels
		<< " (expected, boundary-inclusive: " << square_expected << ")" << endl;
	bool square_area_ok = square_pixels == square_expected;
	cout << "  " << (square_area_ok ? "PASS" : "FAIL") << endl;
	if (!square_area_ok) failures++;

	float l_analytic = shoelaceArea(l_x, l_y);
	int l_pixels = (int)lshape.getOwnedPixels().size();
	//hand-derived: bottom band x in [0,10] y in [0,5] inclusive (11*6=66) plus the new rows
	//the left arm adds above it, x in [0,5] y in [6,10] inclusive (6*5=30, y=5 row already
	//counted in the bottom band) = 96
	int l_expected = 11 * 6 + 6 * 5;
	cout << "  L-shape: analytic area=" << l_analytic << ", pixel count=" << l_pixels
		<< " (expected, boundary-inclusive: " << l_expected << ")" << endl;
	bool l_area_ok = l_pixels == l_expected;
	cout << "  " << (l_area_ok ? "PASS" : "FAIL") << endl;
	if (!l_area_ok) failures++;

	//--- circle ---
	cout << endl << "=== Circle ===" << endl;
	Circle2D circle(50, 50, 20.f);
	bool circle_ok = circle.contains(50, 50)      //center
		&& circle.contains(50, 69)                //just inside the top (radius 20, so y=30..70)
		&& !circle.contains(50, 29)                //just outside the top
		&& !circle.contains(80, 80);               //far outside
	cout << "  " << (circle_ok ? "PASS" : "FAIL") << ": spot checks" << endl;
	if (!circle_ok) failures++;

	float circle_analytic = 3.14159265f * 20.f * 20.f;
	int circle_pixels = (int)circle.getOwnedPixels().size();
	cout << "  analytic area=" << circle_analytic << ", pixel count=" << circle_pixels << endl;
	bool circle_area_ok = fabs(circle_pixels - circle_analytic) < 0.05f * circle_analytic;
	cout << "  " << (circle_area_ok ? "PASS" : "FAIL") << endl;
	if (!circle_area_ok) failures++;

	//--- practical use: masked POI queue ---
	cout << endl << "=== Practical use: masked POI queue from getOwnedPixels() ===" << endl;
	vector<POI2D> poi_queue;
	for (auto& px : circle.getOwnedPixels())
	{
		poi_queue.push_back(POI2D(px));
	}
	bool all_inside = true;
	for (auto& poi : poi_queue)
	{
		if (!circle.contains((int)poi.x, (int)poi.y)) all_inside = false;
	}
	cout << "  built " << poi_queue.size() << " POIs directly from the circle's owned pixels; all inside: " << (all_inside ? "yes" : "no") << endl;
	cout << "  " << (all_inside && poi_queue.size() == circle.getOwnedPixels().size() ? "PASS" : "FAIL") << endl;
	if (!(all_inside && poi_queue.size() == circle.getOwnedPixels().size())) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
