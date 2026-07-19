/*
 Local build-verification test for RegionWithHoles2D, "topology-aware ROI"
 phase 1a (point-membership only -- see oc_shape.h's scope note; connectivity-
 aware subset/strain clipping, phase 1b/2, was scoped separately and is not
 currently planned).

 1. Spot checks: a square with a circular hole in the middle -- points inside
    the hole excluded, points elsewhere in the square included, points outside
    the square excluded regardless of the hole.
 2. Area cross-check: RegionWithHoles2D's getOwnedPixels() count should equal
    the outer shape's own pixel count minus the hole's own pixel count
    exactly, since a hole fully contained within the outer shape can't
    overlap itself. Independent Polygon2D/Circle2D instances are used for the
    "outer area"/"hole area" reference counts vs. the ones actually moved into
    RegionWithHoles2D, since ownership transfers via unique_ptr.
 3. A hole that pokes outside the outer boundary: only the overlapping portion
    should be excluded (contains() is just "outer AND NOT any hole", so this
    should fall out naturally with no special-casing).
 4. Multiple holes in one region.
 5. Practical use: masked POI queue built from getOwnedPixels(), confirming
    every POI is inside the outer shape and outside every hole.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	//--- square (10,10)-(20,20) with a circular hole at its center ---
	cout << "=== Square with a circular hole ===" << endl;
	vector<int> sq_x = { 10, 20, 20, 10 };
	vector<int> sq_y = { 10, 10, 20, 20 };

	auto outer_for_region = std::make_unique<Polygon2D>(sq_x, sq_y);
	vector<std::unique_ptr<Shape2D>> holes_for_region;
	holes_for_region.push_back(std::make_unique<Circle2D>(15, 15, 3.f));
	RegionWithHoles2D square_with_hole(std::move(outer_for_region), std::move(holes_for_region));

	bool spot_ok = !square_with_hole.contains(15, 15)   //hole center: excluded
		&& !square_with_hole.contains(15, 13)           //still inside the hole (radius 3): excluded
		&& square_with_hole.contains(15, 10)            //square edge, well outside the hole: included
		&& square_with_hole.contains(11, 11)            //square corner region: included
		&& !square_with_hole.contains(5, 15)            //outside the square entirely: excluded
		&& !square_with_hole.contains(15, 25);           //outside the square entirely: excluded
	cout << "  " << (spot_ok ? "PASS" : "FAIL") << ": hole excluded, rest of square included, outside-square excluded" << endl;
	if (!spot_ok) failures++;

	//--- area cross-check ---
	cout << endl << "=== Area cross-check (getOwnedPixels count == outer count - hole count) ===" << endl;
	Polygon2D outer_reference(sq_x, sq_y);
	Circle2D hole_reference(15, 15, 3.f);
	int outer_pixels = (int)outer_reference.getOwnedPixels().size();
	int hole_pixels = (int)hole_reference.getOwnedPixels().size();
	int region_pixels = (int)square_with_hole.getOwnedPixels().size();
	cout << "  outer=" << outer_pixels << " hole=" << hole_pixels << " region=" << region_pixels
		<< " (expected " << (outer_pixels - hole_pixels) << ")" << endl;
	bool area_ok = region_pixels == outer_pixels - hole_pixels;
	cout << "  " << (area_ok ? "PASS" : "FAIL") << endl;
	if (!area_ok) failures++;

	//--- hole partially outside the outer boundary ---
	cout << endl << "=== Hole partially outside the outer boundary ===" << endl;
	auto outer2 = std::make_unique<Polygon2D>(sq_x, sq_y);
	vector<std::unique_ptr<Shape2D>> holes2;
	holes2.push_back(std::make_unique<Circle2D>(20, 15, 4.f)); //centered ON the right edge, pokes outside
	RegionWithHoles2D edge_hole_region(std::move(outer2), std::move(holes2));

	bool edge_hole_ok = !edge_hole_region.contains(18, 15)  //inside square AND inside the hole: excluded
		&& edge_hole_region.contains(12, 15)                //inside square, outside the hole: included
		&& !edge_hole_region.contains(23, 15);               //outside the square (hole extends here, but square doesn't): excluded
	cout << "  " << (edge_hole_ok ? "PASS" : "FAIL") << ": only the overlapping portion is excluded" << endl;
	if (!edge_hole_ok) failures++;

	//--- multiple holes ---
	cout << endl << "=== Multiple holes ===" << endl;
	auto outer3 = std::make_unique<Polygon2D>(sq_x, sq_y);
	vector<std::unique_ptr<Shape2D>> holes3;
	holes3.push_back(std::make_unique<Circle2D>(13, 13, 1.f));
	holes3.push_back(std::make_unique<Circle2D>(17, 17, 1.f));
	RegionWithHoles2D multi_hole_region(std::move(outer3), std::move(holes3));

	bool multi_hole_ok = !multi_hole_region.contains(13, 13)  //first hole: excluded
		&& !multi_hole_region.contains(17, 17)                //second hole: excluded
		&& multi_hole_region.contains(15, 15);                //between the two holes: included
	cout << "  " << (multi_hole_ok ? "PASS" : "FAIL") << ": both holes excluded independently" << endl;
	if (!multi_hole_ok) failures++;

	//--- practical use: masked POI queue ---
	cout << endl << "=== Practical use: masked POI queue from getOwnedPixels() ===" << endl;
	vector<POI2D> poi_queue;
	for (auto& px : square_with_hole.getOwnedPixels())
	{
		poi_queue.push_back(POI2D(px));
	}
	bool all_valid = true;
	for (auto& poi : poi_queue)
	{
		if (!square_with_hole.contains((int)poi.x, (int)poi.y)) all_valid = false;
	}
	cout << "  built " << poi_queue.size() << " POIs directly from the region's owned pixels; all valid: " << (all_valid ? "yes" : "no") << endl;
	bool practical_ok = all_valid && poi_queue.size() == square_with_hole.getOwnedPixels().size();
	cout << "  " << (practical_ok ? "PASS" : "FAIL") << endl;
	if (!practical_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
