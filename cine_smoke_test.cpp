/*
 Local build-verification test for Cine2D (.cine file I/O), issue #4.
 Decodes known frames from three real .cine files (8-bit, 10/12-bit-packed,
 16-bit -- the three bit-depth paths hypercine handles) and compares them
 against known-good reference TIFFs bundled with dicengine/hypercine's own
 test suite. Then runs a real FFTCC+ICGN correlation between two decoded
 cine frames to confirm the resulting Image2D is usable by the rest of
 OpenCorr, not just structurally valid.
*/

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

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

	//--- regression: a truncated .cine file must fail cleanly, not silently seek to
	//garbage (hypercine's HyperCine::read_header(), deps/hypercine/hypercine.cpp) ---
	//
	//Truncating a real bundled .cine file lands on an EARLIER, pre-existing sanity check
	//(bitmap_header_.size_image*header_.image_count <= file_size) for essentially any
	//meaningful truncation, since the real files here are dominated by image-data bytes --
	//that check alone already throws for those cases, before this fix's own guards would
	//ever run. To actually exercise this fix (frame_rate and the image-offsets loop, both
	//previously read with no check at all), a minimal SYNTHETIC .cine byte layout is
	//hand-built below with just enough of the header/bitmap-header fields to satisfy every
	//pre-existing check, then truncated exactly at the point where frame_rate would be
	//read -- isolating the specific gap this fix closes.
	cout << endl << "=== Regression: truncated cine file fails cleanly ===" << endl;
	{
		std::vector<char> buf;
		auto push16 = [&](uint16_t v) { buf.insert(buf.end(), (char*)&v, (char*)&v + 2); };
		auto push32 = [&](uint32_t v) { buf.insert(buf.end(), (char*)&v, (char*)&v + 4); };
		auto push64 = [&](uint64_t v) { buf.insert(buf.end(), (char*)&v, (char*)&v + 8); };

		//CINE header (44 bytes total, matching hypercine's own HEADER_SIZE macro)
		push16(0);        //type
		push16(44);       //header_size == test_size (36 fixed fields + 8-byte trigger_time)
		push16(0);        //compression, must be 0
		push16(1);        //version, must be 1
		push32(0);        //first_movie_image
		push32(2);        //total_image_count
		push32(0);        //first_image_no
		push32(2);        //image_count -- >=2 so the unrelated "must have at least two images" check doesn't fire first
		push32(44);       //off_image_header, must == test_size (44)
		push32(84);       //off_setup -- points exactly at EOF of this synthetic file
		push32(84);       //off_image_offsets -- also exactly at EOF
		push64(0);        //trigger_time (TIME64, 8 bytes)

		//BITMAP header (40 bytes, matching header_test_size)
		push32(40);       //size == header_test_size
		push32(4);        //width
		push32(4);        //height
		push16(1);        //planes
		push16(8);        //bit_count, must be 8 or 16
		push32(0);        //compression
		push32(16);       //size_image (4*4*1 byte/px) -- size_image*image_count=32 <= file_size(84)
		push32(0);        //x_pixels_per_meter
		push32(0);        //y_pixels_per_meter
		push32(0);        //clr_used, must be 0
		push32(0);        //clr_important

		bool synth_layout_ok = buf.size() == 84;
		cout << "  synthetic header built: " << buf.size() << " bytes (expected 84)" << endl;

		string truncated_path = "examples/cine/truncated_smoke_test.cine";
		std::ofstream truncated_file(truncated_path, std::ios::binary);
		truncated_file.write(buf.data(), (std::streamsize)buf.size());
		truncated_file.close();

		bool threw_cleanly = false;
		string what_message;
		try
		{
			Cine2D truncated_cine(truncated_path);
		}
		catch (std::exception& e)
		{
			threw_cleanly = true;
			what_message = e.what();
		}
		catch (...)
		{
			//some non-std::exception was thrown -- still "failed cleanly" in the sense of
			//not silently succeeding, but not the descriptive exception the fix produces
		}

		remove(truncated_path.c_str());

		cout << "  threw=" << threw_cleanly << " what=\"" << what_message << "\"" << endl;
		bool truncation_ok = synth_layout_ok && threw_cleanly && what_message.find("truncated cine file") != string::npos;
		cout << "  " << (truncation_ok ? "PASS" : "FAIL")
			<< ": truncated file throws a clear, descriptive exception instead of silently proceeding" << endl;
		if (!truncation_ok) failures++;
	}

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	return failures == 0 ? 0 : 1;
}
