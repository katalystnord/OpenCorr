#include <cstdio>
#include <fstream>
#include <iostream>
#include "opencorr.h"
using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	// --- round-trip test: new format ---
	vector<POI2D> poi_queue;
	POI2D poi(10.f, 20.f);
	poi.deformation.u = 1.5f; poi.deformation.v = -0.5f;
	poi.result.u0 = 1.f; poi.result.v0 = -0.4f; poi.result.zncc = 0.98f;
	poi.result.iteration = 3.f; poi.result.convergence = 0.001f; poi.result.feature = 5.f;
	poi.result.sigma = 0.012f; poi.result.beta = 0.34f;
	poi.strain.exx = 0.01f; poi.strain.eyy = 0.02f; poi.strain.exy = 0.03f;
	poi.subset_radius.x = 16.f; poi.subset_radius.y = 16.f;
	poi_queue.push_back(poi);

	IO2D io;
	io.setDelimiter(",");
	io.setPath("/tmp/oc_io_roundtrip.csv");
	io.saveTable2D(poi_queue);

	vector<POI2D> loaded = io.loadTable2D();
	if (loaded.size() != 1) { cout << "FAIL: wrong row count " << loaded.size() << endl; failures++; }
	else {
		POI2D& p = loaded[0];
		bool ok = fabs(p.result.sigma - 0.012f) < 1e-4f && fabs(p.result.beta - 0.34f) < 1e-4f
			&& fabs(p.strain.exx - 0.01f) < 1e-4f && fabs(p.strain.eyy - 0.02f) < 1e-4f
			&& fabs(p.strain.exy - 0.03f) < 1e-4f
			&& fabs(p.subset_radius.x - 16.f) < 1e-4f && fabs(p.subset_radius.y - 16.f) < 1e-4f;
		cout << "round-trip: sigma=" << p.result.sigma << " beta=" << p.result.beta
			<< " exx=" << p.strain.exx << " eyy=" << p.strain.eyy << " exy=" << p.strain.exy
			<< " subset_rx=" << p.subset_radius.x << " subset_ry=" << p.subset_radius.y << endl;
		cout << (ok ? "PASS" : "FAIL") << ": new-format round-trip" << endl;
		if (!ok) failures++;
	}

	// header check
	ifstream f("/tmp/oc_io_roundtrip.csv");
	string header; getline(f, header);
	cout << "header: " << header << endl;
	bool header_ok = header.find("sigma") != string::npos && header.find("beta") != string::npos
		&& header.find("exx") != string::npos;
	cout << (header_ok ? "PASS" : "FAIL") << ": header contains sigma/beta" << endl;
	if (!header_ok) failures++;

	// --- legacy format test: hand-written 15-column CSV (pre-sigma/beta) ---
	ofstream legacy("/tmp/oc_io_legacy.csv");
	legacy << "x,y,u,v,u0,v0,ZNCC,iteration,convergence,feature,exx,eyy,exy,subset_rx,subset_ry" << endl;
	legacy << "10,20,1.5,-0.5,1,-0.4,0.98,3,0.001,5,0.01,0.02,0.03,16,16" << endl;
	legacy.close();

	IO2D io_legacy;
	io_legacy.setDelimiter(",");
	io_legacy.setPath("/tmp/oc_io_legacy.csv");
	vector<POI2D> legacy_loaded = io_legacy.loadTable2D();
	if (legacy_loaded.size() != 1) { cout << "FAIL: legacy wrong row count " << legacy_loaded.size() << endl; failures++; }
	else {
		POI2D& p = legacy_loaded[0];
		bool ok = p.result.sigma == -1.f && p.result.beta == 0.f
			&& fabs(p.strain.exx - 0.01f) < 1e-4f && fabs(p.strain.eyy - 0.02f) < 1e-4f
			&& fabs(p.strain.exy - 0.03f) < 1e-4f
			&& fabs(p.subset_radius.x - 16.f) < 1e-4f && fabs(p.subset_radius.y - 16.f) < 1e-4f;
		cout << "legacy load: sigma=" << p.result.sigma << " beta=" << p.result.beta
			<< " exx=" << p.strain.exx << " eyy=" << p.strain.eyy << " exy=" << p.strain.exy
			<< " subset_rx=" << p.subset_radius.x << " subset_ry=" << p.subset_radius.y << endl;
		cout << (ok ? "PASS" : "FAIL") << ": legacy 15-column CSV loads correctly, sigma/beta sentinels set" << endl;
		if (!ok) failures++;
	}

	// --- loadDeformationTable2D(): previously had no counterpart to saveDeformationTable2D(),
	// so the full 12-component deformation vector could be saved but never loaded back ---
	POI2D poi_def(30.f, 40.f);
	poi_def.deformation.u = 1.1f; poi_def.deformation.ux = 0.01f; poi_def.deformation.uy = 0.02f;
	poi_def.deformation.uxx = 0.001f; poi_def.deformation.uxy = 0.002f; poi_def.deformation.uyy = 0.003f;
	poi_def.deformation.v = -0.9f; poi_def.deformation.vx = 0.04f; poi_def.deformation.vy = 0.05f;
	poi_def.deformation.vxx = 0.004f; poi_def.deformation.vxy = 0.005f; poi_def.deformation.vyy = 0.006f;
	poi_def.subset_radius.x = 12.f; poi_def.subset_radius.y = 14.f;
	vector<POI2D> poi_def_queue = { poi_def };

	IO2D io_def;
	io_def.setDelimiter(",");
	io_def.setPath("/tmp/oc_io_deformation.csv");
	io_def.saveDeformationTable2D(poi_def_queue);
	vector<POI2D> loaded_def = io_def.loadDeformationTable2D();
	if (loaded_def.size() != 1) { cout << "FAIL: loadDeformationTable2D wrong row count " << loaded_def.size() << endl; failures++; }
	else {
		POI2D& p = loaded_def[0];
		bool def_ok = true;
		for (int i = 0; i < 12; i++)
		{
			if (fabs(p.deformation.p[i] - poi_def.deformation.p[i]) > 1e-4f) def_ok = false;
		}
		def_ok = def_ok && fabs(p.subset_radius.x - 12.f) < 1e-4f && fabs(p.subset_radius.y - 14.f) < 1e-4f;
		cout << (def_ok ? "PASS" : "FAIL") << ": loadDeformationTable2D() round-trips all 12 deformation components" << endl;
		if (!def_ok) failures++;
	}

	// --- malformed CSV row must be skipped, not crash the whole load ---
	ofstream malformed("/tmp/oc_io_malformed.csv");
	malformed << "x,y,u,v,u0,v0,ZNCC,iteration,convergence,feature,sigma,beta,exx,eyy,exy,subset_rx,subset_ry" << endl;
	malformed << "10,20,1.5,-0.5,1,-0.4,0.98,3,0.001,5,0.012,0.34,0.01,0.02,0.03,16,16" << endl;
	malformed << "" << endl; //blank line
	malformed << "5,6" << endl; //truncated row, far too few fields
	malformed << "10,20,1.5,-0.5,1,-0.4,0.98,3,0.001,5,0.012,0.34,0.01,0.02,0.03,16,16" << endl; //well-formed again
	malformed.close();

	IO2D io_malformed;
	io_malformed.setDelimiter(",");
	io_malformed.setPath("/tmp/oc_io_malformed.csv");
	vector<POI2D> malformed_loaded = io_malformed.loadTable2D(); //must not crash
	bool malformed_ok = malformed_loaded.size() == 2; //the two well-formed rows, malformed ones skipped
	cout << "loaded " << malformed_loaded.size() << " POIs from a file with a blank and a truncated row mixed in" << endl;
	cout << (malformed_ok ? "PASS" : "FAIL") << ": malformed rows are skipped, not a crash or silent out-of-bounds read" << endl;
	if (!malformed_ok) failures++;

	// --- saveMap2D: a failed POI's sentinel must not leak into the raster as a real value ---
	vector<POI2D> map_queue;
	POI2D good_poi(1.f, 1.f);
	good_poi.result.zncc = 0.95f;
	good_poi.deformation.u = 2.5f;
	map_queue.push_back(good_poi);
	POI2D failed_poi(2.f, 1.f);
	failed_poi.result.zncc = (float)STATUS_HESSIAN_SINGULAR; //a real solver failure sentinel
	failed_poi.deformation.u = 999.f; //stale/garbage, as an unsolved POI's would be
	map_queue.push_back(failed_poi);

	IO2D io_map;
	io_map.setDelimiter(",");
	io_map.setHeight(4);
	io_map.setWidth(4);
	io_map.setPath("/tmp/oc_io_map.csv");
	io_map.saveMap2D(map_queue, u);

	ifstream map_file("/tmp/oc_io_map.csv");
	vector<string> map_lines;
	string map_line;
	while (getline(map_file, map_line)) map_lines.push_back(map_line);
	//row 1 (y=1) should hold "2.5" at column x=1 (the good POI) and "nan" at column x=2 (the failed one)
	bool map_ok = map_lines.size() > 1 && map_lines[1].find("2.5") != string::npos
		&& (map_lines[1].find("nan") != string::npos || map_lines[1].find("-nan") != string::npos);
	cout << "row for y=1: " << (map_lines.size() > 1 ? map_lines[1] : "<missing>") << endl;
	cout << (map_ok ? "PASS" : "FAIL") << ": failed POI's sentinel becomes NaN in the map, not a raw -8 next to real values" << endl;
	if (!map_ok) failures++;

	// --- saveMap2D/saveMap3D: displacement-gradient OutputVariable codes (u_x etc.) must
	// actually be written, not silently no-op through the switch's default case ---
	good_poi.deformation.ux = 3.25f; //nonzero and distinct from the map's own 0.0 default,
	map_queue[0] = good_poi;         //so "wrote the real value" is distinguishable from "no-op"
	remove("/tmp/oc_io_map.csv"); //a stale file from an earlier run must not fake a PASS
	io_map.saveMap2D(map_queue, u_x);
	ifstream map_ux_file("/tmp/oc_io_map.csv");
	vector<string> map_ux_lines;
	string map_ux_line;
	while (getline(map_ux_file, map_ux_line)) map_ux_lines.push_back(map_ux_line);
	bool map_2d_ux_ok = map_ux_lines.size() > 1 && map_ux_lines[1].find("3.25") != string::npos;
	cout << (map_2d_ux_ok ? "PASS" : "FAIL") << ": saveMap2D(u_x) writes the displacement gradient, not a silent no-op" << endl;
	if (!map_2d_ux_ok) failures++;

	POI3D map3d_poi(1.f, 1.f, 1.f);
	map3d_poi.result.zncc = 0.9f;
	map3d_poi.deformation.ux = 4.75f;
	vector<POI3D> map3d_queue = { map3d_poi };
	IO3D io_map3d;
	io_map3d.setDelimiter(",");
	io_map3d.setDimX(4); io_map3d.setDimY(4); io_map3d.setDimZ(4);
	io_map3d.setPath("/tmp/oc_io_map3d.csv");
	remove("/tmp/oc_io_map3d.csv"); //a stale file from an earlier run must not fake a PASS
	io_map3d.saveMap3D(map3d_queue, u_x);
	ifstream map3d_file("/tmp/oc_io_map3d.csv");
	string map3d_contents((istreambuf_iterator<char>(map3d_file)), istreambuf_iterator<char>());
	bool map_3d_ux_ok = map3d_contents.find("4.75") != string::npos;
	cout << (map_3d_ux_ok ? "PASS" : "FAIL") << ": saveMap3D(u_x) writes the displacement gradient, not a silent no-op" << endl;
	if (!map_3d_ux_ok) failures++;

	// --- saveMatrixBin/loadMatrixBin: binary mode must match, and a malformed/truncated
	// file must return an empty queue instead of crashing ---
	POI3D poi3d(1.f, 2.f, 3.f);
	poi3d.deformation.u = 0.5f; poi3d.deformation.v = -0.3f; poi3d.deformation.w = 0.1f;
	poi3d.result.zncc = 0.9f; poi3d.result.convergence = 0.002f;
	vector<POI3D> poi3d_queue = { poi3d };

	IO3D io3d;
	io3d.setDelimiter(",");
	io3d.setPath("/tmp/oc_io_matrix.bin");
	io3d.setDimX(10); io3d.setDimY(10); io3d.setDimZ(10);
	io3d.saveMatrixBin(poi3d_queue);

	IO3D io3d_load;
	io3d_load.setPath("/tmp/oc_io_matrix.bin");
	vector<POI3D> loaded3d = io3d_load.loadMatrixBin();
	bool bin_ok = loaded3d.size() == 1
		&& fabs(loaded3d[0].x - 1.f) < 1e-4f && fabs(loaded3d[0].deformation.u - 0.5f) < 1e-4f
		&& fabs(loaded3d[0].result.zncc - 0.9f) < 1e-4f;
	cout << (bin_ok ? "PASS" : "FAIL") << ": saveMatrixBin/loadMatrixBin binary round-trip" << endl;
	if (!bin_ok) failures++;

	//a nonexistent file must return an empty queue, not fall through to read() on a closed stream
	IO3D io3d_missing;
	io3d_missing.setPath("/tmp/oc_io_matrix_does_not_exist.bin");
	vector<POI3D> missing_result = io3d_missing.loadMatrixBin();
	cout << (missing_result.empty() ? "PASS" : "FAIL") << ": loadMatrixBin() on a missing file returns an empty queue, not a crash" << endl;
	if (!missing_result.empty()) failures++;

	//a truncated file (valid header claiming 1 POI, but the data payload cut short) must
	//return empty rather than reading garbage/uninitialized memory into the result
	ifstream src("/tmp/oc_io_matrix.bin", ios::binary);
	vector<char> full_bytes((istreambuf_iterator<char>(src)), istreambuf_iterator<char>());
	src.close();
	ofstream truncated("/tmp/oc_io_matrix_truncated.bin", ios::binary);
	truncated.write(full_bytes.data(), (long)full_bytes.size() - 8); //chop off the last two floats
	truncated.close();

	IO3D io3d_truncated;
	io3d_truncated.setPath("/tmp/oc_io_matrix_truncated.bin");
	vector<POI3D> truncated_result = io3d_truncated.loadMatrixBin();
	cout << (truncated_result.empty() ? "PASS" : "FAIL") << ": loadMatrixBin() on a truncated file returns an empty queue, not partially-garbage data" << endl;
	if (!truncated_result.empty()) failures++;

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << ")" << endl;
	return failures == 0 ? 0 : 1;
}
