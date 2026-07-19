/*
 Local build-verification test for the StatusFlag taxonomy + Hessian-
 singularity guard, issue #13.

 Three checks:
 1. statusDescription() returns a non-empty, distinct string for every known
    code, and reports success for a non-negative ZNCC.
 2. Motivating case for the new STATUS_HESSIAN_SINGULAR guard: a subset with
    real, finite intensity gradient that varies in x only (uniform along y --
    a stand-in for a real unidirectional-texture surface, e.g. a machined or
    rolled surface with directional lay) makes ICGN2D1's Hessian exactly
    rank-deficient in its y-shape-function block. Before this guard, that
    silently flowed into hessian.inverse() with no invertibility check.
 3. Regression: existing named codes (out-of-bounds subset, NaN incoming
    guess) still fire correctly after the refactor from magic numbers to
    StatusFlag constants, and a normal well-textured subset still succeeds.
*/

#include <iostream>
#include <set>
#include <limits>
#include <cmath>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	//--- check 1: statusDescription() ---
	cout << "=== statusDescription() sanity ===" << endl;
	vector<float> codes = {
		1.f, //success
		(float)STATUS_INSUFFICIENT_FEATURES,
		(float)STATUS_DEGENERATE_INPUT,
		(float)STATUS_INVALID_SUBSET_OR_GUESS,
		(float)STATUS_MAX_ITERATIONS_REACHED,
		(float)STATUS_NAN_IN_RESULT,
		(float)STATUS_RELIABILITY_GUIDED_REJECTED,
		(float)STATUS_SEQUENCE_JUMP_REJECTED,
		(float)STATUS_HESSIAN_SINGULAR,
	};
	std::set<std::string> seen_descriptions;
	bool descriptions_ok = true;
	for (float code : codes)
	{
		std::string desc = statusDescription(code);
		cout << "  " << code << " -> \"" << desc << "\"" << endl;
		if (desc.empty() || seen_descriptions.count(desc) != 0)
		{
			descriptions_ok = false;
		}
		seen_descriptions.insert(desc);
	}
	if (statusDescription(1.f).find("succeeded") == std::string::npos)
	{
		descriptions_ok = false;
	}
	cout << "  " << (descriptions_ok ? "PASS" : "FAIL") << ": every code has a distinct, non-empty description" << endl;
	if (!descriptions_ok) failures++;

	//isFailureStatus() must agree with statusDescription() on every named code (both are
	//driven by the same table in oc_dic.cpp) -- this is the actual regression test for the
	//bug that motivated moving it out of ReliabilityGuided2D's own local
	//isSolverFailureSentinel(), which had silently gone stale missing
	//STATUS_RELIABILITY_GUIDED_REJECTED and STATUS_SEQUENCE_JUMP_REJECTED
	bool sentinel_ok = true;
	for (float code : codes)
	{
		if (code >= 0.f) continue; //the leading 1.f success placeholder
		if (!isFailureStatus(code))
		{
			cout << "  isFailureStatus(" << code << ") incorrectly returned false" << endl;
			sentinel_ok = false;
		}
	}
	if (isFailureStatus(1.f) || isFailureStatus(0.f) || isFailureStatus(-0.3f))
	{
		cout << "  isFailureStatus() incorrectly flagged a genuine (non-sentinel) value as a failure" << endl;
		sentinel_ok = false;
	}
	cout << "  " << (sentinel_ok ? "PASS" : "FAIL") << ": isFailureStatus() recognizes every named code, and only those" << endl;
	if (!sentinel_ok) failures++;

	//--- check 2: Hessian-singularity guard ---
	cout << endl << "=== Motivating case: unidirectional-texture subset triggers STATUS_HESSIAN_SINGULAR ===" << endl;
	int img_size = 80;
	Image2D dir_ref(img_size, img_size);
	Image2D dir_tar(img_size, img_size);

	//intensity varies only along x (a sinusoid); every row is identical, so the
	//finite-difference gradient along y is exactly zero everywhere -- the
	//uy/vy-related rows and columns of the 6x6 Hessian are then exactly zero,
	//making it exactly rank-deficient (singular), not merely ill-conditioned
	for (int y = 0; y < img_size; y++)
	{
		for (int x = 0; x < img_size; x++)
		{
			float val = 128.f + 100.f * sin((float)x * 0.3f);
			uchar pixel = (uchar)max(0.f, min(255.f, val));
			dir_ref.cv_mat.at<uchar>(y, x) = pixel;
			dir_tar.cv_mat.at<uchar>(y, x) = pixel;
		}
	}
	cv::cv2eigen(dir_ref.cv_mat, dir_ref.eg_mat);
	cv::cv2eigen(dir_tar.cv_mat, dir_tar.eg_mat);

	POI2D poi_singular(Point2D(40.f, 40.f));
	ICGN2D1 icgn_singular(20, 20, 0.001f, 20, 1);
	icgn_singular.setImages(dir_ref, dir_tar);
	icgn_singular.prepare();
	icgn_singular.compute(&poi_singular);

	cout << "  zncc=" << poi_singular.result.zncc << " (\"" << statusDescription(poi_singular.result.zncc) << "\")"
		<< " u=" << poi_singular.deformation.u << " v=" << poi_singular.deformation.v << endl;
	bool singular_ok = poi_singular.result.zncc == (float)STATUS_HESSIAN_SINGULAR
		&& !std::isnan(poi_singular.deformation.u) && !std::isnan(poi_singular.deformation.v);
	cout << "  " << (singular_ok ? "PASS" : "FAIL") << ": singular Hessian rejected cleanly instead of silently inverted" << endl;
	if (!singular_ok) failures++;

	//--- check 3: regression -- existing named codes still fire correctly ---
	cout << endl << "=== Regression: existing sentinel behaviors after the magic-number refactor ===" << endl;
	Image2D ref_img("examples/2d_dic/oht_cfrp_0.bmp");
	Image2D tar_img("examples/2d_dic/oht_cfrp_4.bmp");

	ICGN2D1 icgn(16, 16, 0.001f, 20, 1);
	icgn.setImages(ref_img, tar_img);
	icgn.prepare();

	//out-of-bounds subset (POI too close to the image edge)
	POI2D poi_oob(Point2D(2.f, 2.f));
	icgn.compute(&poi_oob);
	bool oob_ok = poi_oob.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  out-of-bounds subset -> zncc=" << poi_oob.result.zncc
		<< " (\"" << statusDescription(poi_oob.result.zncc) << "\") " << (oob_ok ? "PASS" : "FAIL") << endl;
	if (!oob_ok) failures++;

	//NaN incoming initial guess
	POI2D poi_nan_guess(Point2D(100.f, 100.f));
	poi_nan_guess.deformation.u = std::numeric_limits<float>::quiet_NaN();
	icgn.compute(&poi_nan_guess);
	bool nan_guess_ok = poi_nan_guess.result.zncc == (float)STATUS_INVALID_SUBSET_OR_GUESS;
	cout << "  NaN incoming guess -> zncc=" << poi_nan_guess.result.zncc
		<< " (\"" << statusDescription(poi_nan_guess.result.zncc) << "\") " << (nan_guess_ok ? "PASS" : "FAIL") << endl;
	if (!nan_guess_ok) failures++;

	//normal well-textured subset still succeeds
	POI2D poi_ok(Point2D(100.f, 100.f));
	icgn.compute(&poi_ok);
	bool normal_ok = poi_ok.result.zncc > 0.9f;
	cout << "  well-textured subset -> zncc=" << poi_ok.result.zncc
		<< " (\"" << statusDescription(poi_ok.result.zncc) << "\") " << (normal_ok ? "PASS" : "FAIL") << endl;
	if (!normal_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
