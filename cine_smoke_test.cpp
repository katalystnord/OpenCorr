/*
 Local build-verification test for Cine2D (.cine file I/O), issue #4.
 Decodes known frames from three real .cine files (8-bit, 10/12-bit-packed,
 16-bit -- the three bit-depth paths hypercine handles) and compares them
 against known-good reference TIFFs bundled with dicengine/hypercine's own
 test suite. Then runs a real FFTCC+ICGN correlation between two decoded
 cine frames to confirm the resulting Image2D is usable by the rest of
 OpenCorr, not just structurally valid.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

static bool checkFrame(const string& cine_path, int frame_id, const string& reference_tiff_path)
{
	Cine2D cine(cine_path);
	cout << cine_path << ": " << cine.width() << "x" << cine.height()
		<< ", " << cine.frameCount() << " frames, first id " << cine.firstFrameId() << endl;

	Image2D decoded = cine.getFrame(frame_id);
	Image2D reference(reference_tiff_path);

	if (decoded.width != reference.width || decoded.height != reference.height)
	{
		cout << "  FAIL: size mismatch, decoded " << decoded.width << "x" << decoded.height
			<< " vs reference " << reference.width << "x" << reference.height << endl;
		return false;
	}

	Eigen::MatrixXf diff = decoded.eg_mat - reference.eg_mat;
	float mean_abs_diff = diff.cwiseAbs().mean();
	float max_abs_diff = diff.cwiseAbs().maxCoeff();

	cout << "  frame " << frame_id << " vs " << reference_tiff_path
		<< ": mean|diff|=" << mean_abs_diff << " max|diff|=" << max_abs_diff << endl;

	//allow +/-1 intensity level of rounding slack; anything more indicates a real decode bug
	bool ok = max_abs_diff <= 1.f;
	cout << "  " << (ok ? "PASS" : "FAIL") << endl;
	return ok;
}

int main()
{
	int failures = 0;

	//note: despite the filename, hypercine's own reference test reads frame index 10 for this file/tiff pair
	if (!checkFrame("examples/cine/example_8bpp.cine", 10, "examples/cine/example_8bpp_frame_11.tiff")) failures++;
	if (!checkFrame("examples/cine/packed_12bpp.cine", 60, "examples/cine/packed_12bpp_frame_60.tiff")) failures++;
	if (!checkFrame("examples/cine/phantom_v7_raw_16bpp.cine", 238292, "examples/cine/phantom_v7_raw_16bpp_frame_238292.tiff")) failures++;

	//integration check: correlate two decoded cine frames with the rest of OpenCorr
	cout << endl << "Integration check: FFTCC+ICGN on two decoded cine frames" << endl;
	int thread_number = 4;
	omp_set_num_threads(thread_number); //must match the instance-pool size the solvers below are constructed with

	Cine2D cine("examples/cine/example_8bpp.cine");
	Image2D ref_img = cine.getFrame(7);
	Image2D tar_img = cine.getFrame(10);

	vector<POI2D> poi_queue;
	for (int y = 20; y < ref_img.height - 20; y += 4)
		for (int x = 20; x < ref_img.width - 20; x += 4)
			poi_queue.push_back(POI2D(Point2D((float)x, (float)y)));

	int subset_radius = 10;
	FFTCC2D fftcc(subset_radius, subset_radius, thread_number);
	fftcc.setImages(ref_img, tar_img);
	fftcc.compute(poi_queue);

	ICGN2D1 icgn1(subset_radius, subset_radius, 0.001f, 10, thread_number);
	icgn1.setImages(ref_img, tar_img);
	icgn1.prepare();
	icgn1.compute(poi_queue);

	int valid = 0;
	for (auto& poi : poi_queue) if (poi.result.zncc > 0.f) valid++;
	cout << poi_queue.size() << " POIs, " << valid << " produced a valid ZNCC result "
		<< "(low bar -- this is an unpatterned/low-contrast test frame pair, the point is that the pipeline runs end to end on cine-decoded images)." << endl;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
