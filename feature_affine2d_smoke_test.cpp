/*
 Local build-verification test for FeatureAffine2D, added alongside the
 FeatureAffine2D/3D deduplication refactor (issue #17, Tier 3) -- prior to
 this test, FeatureAffine2D had no dedicated smoke coverage at all (only
 FeatureAffine3D did), even though the RANSAC/affine-fit core the two share
 is exactly the code most at risk from that refactor.

 Ground truth: a known full 2D affine transform with distinct, nonzero
 ux/uy/vx/vy, applied to a synthetic keypoint cloud with no noise, so
 recovery should be exact (to floating-point precision).
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	cout << "=== FeatureAffine2D: recover a known full 2D affine transform ===" << endl;

	//true transform: distinct, nonzero ux/uy/vx/vy/u/v so any axis mixup would be caught
	float true_u = 1.5f, true_ux = 0.2f, true_uy = 0.1f;
	float true_v = -2.0f, true_vx = 0.05f, true_vy = 0.3f;

	vector<Point2D> ref_kp, tar_kp;
	for (float x = -5.f; x <= 5.f; x += 1.f)
	{
		for (float y = -5.f; y <= 5.f; y += 1.f)
		{
			ref_kp.push_back(Point2D(x, y));
			Point2D tar_pt;
			tar_pt.x = (1.f + true_ux) * x + true_uy * y + true_u;
			tar_pt.y = true_vx * x + (1.f + true_vy) * y + true_v;
			tar_kp.push_back(tar_pt);
		}
	}

	FeatureAffine2D feature_affine(10, 10, 2);
	feature_affine.setSearch(15.f, 7);
	RansacConfig ransac_config;
	ransac_config.trial_number = 10;
	ransac_config.sample_mumber = 3;
	ransac_config.error_threshold = 1.f;
	feature_affine.setRansacConfig(ransac_config);
	feature_affine.setKeypointPair(ref_kp, tar_kp);
	feature_affine.prepare();

	POI2D poi(Point2D(0.f, 0.f));
	feature_affine.compute(&poi);

	cout << "  recovered: u=" << poi.deformation.u << " ux=" << poi.deformation.ux << " uy=" << poi.deformation.uy
		<< " v=" << poi.deformation.v << " vx=" << poi.deformation.vx << " vy=" << poi.deformation.vy << endl;
	cout << "  ground truth: u=" << true_u << " ux=" << true_ux << " uy=" << true_uy
		<< " v=" << true_v << " vx=" << true_vx << " vy=" << true_vy << endl;

	bool ok = fabs(poi.deformation.u - true_u) < 1e-3f
		&& fabs(poi.deformation.ux - true_ux) < 1e-3f
		&& fabs(poi.deformation.uy - true_uy) < 1e-3f
		&& fabs(poi.deformation.v - true_v) < 1e-3f
		&& fabs(poi.deformation.vx - true_vx) < 1e-3f
		&& fabs(poi.deformation.vy - true_vy) < 1e-3f
		&& poi.result.zncc == 0.f;
	cout << "  " << (ok ? "PASS" : "FAIL") << ": all six 1st-order deformation terms recovered correctly" << endl;
	if (!ok) failures++;

	//--- regression: self-adaptive subset path (the 2D-only branch, not shared with 3D) ---
	cout << endl << "=== FeatureAffine2D: self-adaptive subset sizing ===" << endl;
	//self-adaptive mode reads ref_img->width/height to clamp the estimated subset
	//placement -- a real Image2D must be attached via setImages() even though this
	//module otherwise operates purely on the keypoint clouds, not pixel content
	Image2D dummy_ref(64, 64), dummy_tar(64, 64);
	FeatureAffine2D feature_affine_adaptive(10, 10, 2);
	feature_affine_adaptive.setImages(dummy_ref, dummy_tar);
	feature_affine_adaptive.setSelfAdaptive(true);
	feature_affine_adaptive.setSearch(15.f, 7);
	feature_affine_adaptive.setSubsetAdjustment(14, 3);
	feature_affine_adaptive.setRansacConfig(ransac_config);
	feature_affine_adaptive.setKeypointPair(ref_kp, tar_kp);
	feature_affine_adaptive.prepare();

	POI2D poi_adaptive(Point2D(0.f, 0.f));
	feature_affine_adaptive.compute(&poi_adaptive);

	bool adaptive_ok = fabs(poi_adaptive.deformation.u - true_u) < 1e-3f
		&& fabs(poi_adaptive.deformation.vy - true_vy) < 1e-3f
		&& poi_adaptive.subset_radius.x >= 3 && poi_adaptive.subset_radius.y >= 3;
	cout << "  recovered: u=" << poi_adaptive.deformation.u << " vy=" << poi_adaptive.deformation.vy
		<< " subset_radius=(" << poi_adaptive.subset_radius.x << "," << poi_adaptive.subset_radius.y << ")" << endl;
	cout << "  " << (adaptive_ok ? "PASS" : "FAIL") << ": self-adaptive subset path still recovers the correct transform" << endl;
	if (!adaptive_ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
