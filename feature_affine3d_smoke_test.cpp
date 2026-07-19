/*
 Local build-verification test for FeatureAffine3D, surfaced by a full-branch
 audit: poi->deformation.wy was written twice (a copy-paste slip), so wz was
 never actually assigned -- every 3D/DVC feature-guided initial guess silently
 dropped its z-direction displacement gradient.

 Ground truth: a known full 3D affine transform with distinct, nonzero wx/wy/wz
 (so the bug -- wy colliding with wz's true value, wz staying at zero -- would
 be unambiguous), applied to a synthetic keypoint cloud with no noise, so
 recovery should be exact (to floating-point precision).
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	cout << "=== FeatureAffine3D: recover a known full 3D affine transform ===" << endl;

	//true transform: identity in x/y, but z depends on all three axes --
	//wx=0.2, wy=0.1, wz=0.3, w=2.0, all distinct so the wy/wz mixup would be caught
	float true_wx = 0.2f, true_wy = 0.1f, true_wz = 0.3f, true_w = 2.0f;

	vector<Point3D> ref_kp, tar_kp;
	for (float x = -5.f; x <= 5.f; x += 2.f)
	{
		for (float y = -5.f; y <= 5.f; y += 2.f)
		{
			for (float z = -5.f; z <= 5.f; z += 2.f)
			{
				ref_kp.push_back(Point3D(x, y, z));
				Point3D tar_pt;
				tar_pt.x = x; //identity
				tar_pt.y = y; //identity
				tar_pt.z = true_wx * x + true_wy * y + (1.f + true_wz) * z + true_w;
				tar_kp.push_back(tar_pt);
			}
		}
	}

	FeatureAffine3D feature_affine(10, 10, 10, 2);
	feature_affine.setSearch(15.f, 8);
	RansacConfig ransac_config;
	ransac_config.trial_number = 10;
	ransac_config.sample_mumber = 4;
	ransac_config.error_threshold = 1.f;
	feature_affine.setRansacConfig(ransac_config);
	feature_affine.setKeypointPair(ref_kp, tar_kp);
	feature_affine.prepare();

	POI3D poi(Point3D(0.f, 0.f, 0.f));
	feature_affine.compute(&poi);

	cout << "  recovered: wx=" << poi.deformation.wx << " wy=" << poi.deformation.wy
		<< " wz=" << poi.deformation.wz << " w=" << poi.deformation.w << endl;
	cout << "  ground truth: wx=" << true_wx << " wy=" << true_wy << " wz=" << true_wz << " w=" << true_w << endl;

	bool ok = fabs(poi.deformation.wx - true_wx) < 1e-3f
		&& fabs(poi.deformation.wy - true_wy) < 1e-3f
		&& fabs(poi.deformation.wz - true_wz) < 1e-3f
		&& fabs(poi.deformation.w - true_w) < 1e-3f;
	cout << "  " << (ok ? "PASS" : "FAIL") << ": wy and wz recovered independently and correctly "
		"(pre-fix, wy would read " << (1.f + true_wz) << " and wz would read 0)" << endl;
	if (!ok) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
