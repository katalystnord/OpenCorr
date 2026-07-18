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

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << ")" << endl;
	return failures == 0 ? 0 : 1;
}
