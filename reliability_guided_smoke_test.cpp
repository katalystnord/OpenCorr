/*
 Local build-verification test for ReliabilityGuided2D, issue #7.

 Ground truth: solve a full POI grid independently via FFTCC2D+ICGN2D1 (the
 existing, already-verified path). Then solve the SAME grid a second way:
 FFTCC2D+ICGN2D1 on a single seed POI only, then ReliabilityGuided2D flood-
 fill for everything else (ICGN2D1 only, no FFTCC2D call per POI). If the
 propagation mechanism is correct, the two should converge to matching
 results despite the second pass never running a coarse search on most POIs.
*/

#include <iostream>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

namespace
{
	//deterministic stand-in for a real solver, used below to force the exact race the
	//"visited claimed before solving" bug needed: succeeds (converging to the true,
	//exactly-linear field u_true(x) = 0.5*x, and handing forward a correspondingly
	//correct gradient) only if the incoming Taylor-extrapolated guess is already close
	//to truth, and otherwise reports a low ZNCC without moving from the guess -- a
	//reasonable stand-in for "a bad initial guess doesn't converge."
	class LinearFieldStubSolver : public DIC
	{
	public:
		void prepare() override {}

		void compute(POI2D* poi) override
		{
			float true_u = 0.5f * poi->x;
			float error = std::fabs(poi->deformation.u - true_u);
			if (error < 0.01f)
			{
				poi->deformation.u = true_u;
				poi->deformation.ux = 0.5f; //hands the correct gradient forward on success
				poi->result.zncc = 0.9f;
			}
			else
			{
				//didn't converge -- left at the (wrong) guess, low quality
				poi->result.zncc = 0.3f;
			}
		}

		void compute(std::vector<POI2D>& poi_queue) override
		{
			for (auto& poi : poi_queue) compute(&poi);
		}
	};
}

int main()
{
	int failures = 0;
	int subset_radius = 16;
	int thread_number = 4;
	omp_set_num_threads(thread_number);

	Image2D ref_img("examples/2d_dic/oht_cfrp_0.bmp");
	Image2D tar_img("examples/2d_dic/oht_cfrp_4.bmp");

	int poi_number_x = 30, poi_number_y = 30, grid_space = 3;
	Point2D upper_left(60.f, 60.f);

	auto buildGrid = [&]()
	{
		vector<POI2D> grid;
		for (int r = 0; r < poi_number_y; r++)
			for (int c = 0; c < poi_number_x; c++)
				grid.push_back(POI2D(upper_left + Point2D((float)(c * grid_space), (float)(r * grid_space))));
		return grid;
	};

	//--- ground truth: full grid via FFTCC2D + ICGN2D1, independently, for every POI ---
	cout << "=== Ground truth: full-grid FFTCC2D + ICGN2D1 ===" << endl;
	vector<POI2D> ground_truth = buildGrid();

	FFTCC2D fftcc_full(subset_radius, subset_radius, thread_number);
	fftcc_full.setImages(ref_img, tar_img);
	fftcc_full.compute(ground_truth);

	ICGN2D1 icgn_full(subset_radius, subset_radius, 0.001f, 20, thread_number);
	icgn_full.setImages(ref_img, tar_img);
	icgn_full.prepare();
	icgn_full.compute(ground_truth);

	int converged_full = 0;
	for (auto& poi : ground_truth) if (poi.result.zncc > 0.9f) converged_full++;
	cout << "  " << converged_full << " / " << ground_truth.size() << " POIs converged (ZNCC > 0.9)" << endl;

	//--- second pass: single seed + ReliabilityGuided2D flood-fill ---
	cout << endl << "=== Single-seed flood-fill via ReliabilityGuided2D ===" << endl;
	vector<POI2D> flood = buildGrid();

	int seed_row = poi_number_y / 2, seed_col = poi_number_x / 2;
	int seed_idx = seed_row * poi_number_x + seed_col;
	vector<POI2D> seed_only = { flood[seed_idx] };

	FFTCC2D fftcc_seed(subset_radius, subset_radius, 1);
	fftcc_seed.setImages(ref_img, tar_img);
	fftcc_seed.compute(seed_only);

	ICGN2D1 icgn_seed(subset_radius, subset_radius, 0.001f, 20, 1);
	icgn_seed.setImages(ref_img, tar_img);
	icgn_seed.prepare();
	icgn_seed.compute(&seed_only[0]);
	flood[seed_idx] = seed_only[0];

	cout << "  seed POI (" << flood[seed_idx].x << ", " << flood[seed_idx].y << "): zncc=" << flood[seed_idx].result.zncc << endl;

	//the propagation solver: a SEPARATE single-threaded ICGN2D1 instance (its own
	//prepare()'d gradient/interpolation state), called once per POI by ReliabilityGuided2D
	ICGN2D1 icgn_propagate(subset_radius, subset_radius, 0.001f, 20, 1);
	icgn_propagate.setImages(ref_img, tar_img);
	icgn_propagate.prepare();

	ReliabilityGuided2D rg(poi_number_x, poi_number_y);
	rg.setZnccThreshold(0.8f);
	rg.setDeltaDispTolerance(2.f);
	int accepted = rg.compute(flood, { seed_idx }, icgn_propagate);

	int converged_flood = 0;
	for (auto& poi : flood) if (poi.result.zncc > 0.9f) converged_flood++;
	cout << "  flood-fill accepted " << accepted << " POIs beyond the seed, "
		<< converged_flood << " / " << flood.size() << " total converged (ZNCC > 0.9)" << endl;
	cout << "  (note: FFTCC2D/coarse search was called on exactly 1 POI, not " << flood.size() << ")" << endl;

	//--- compare the two solutions where both converged ---
	cout << endl << "=== Comparison ===" << endl;
	float max_diff_u = 0.f, max_diff_v = 0.f;
	int compared = 0;
	for (size_t i = 0; i < ground_truth.size(); i++)
	{
		if (ground_truth[i].result.zncc <= 0.9f || flood[i].result.zncc <= 0.9f) continue;
		float du = fabs(ground_truth[i].deformation.u - flood[i].deformation.u);
		float dv = fabs(ground_truth[i].deformation.v - flood[i].deformation.v);
		max_diff_u = max(max_diff_u, du);
		max_diff_v = max(max_diff_v, dv);
		compared++;
	}
	cout << "  compared " << compared << " POIs where both converged, max |du|=" << max_diff_u << " max |dv|=" << max_diff_v << endl;

	bool coverage_ok = converged_flood > (int)(0.8f * converged_full);
	bool agreement_ok = compared > 0 && max_diff_u < 0.02f && max_diff_v < 0.02f;
	cout << "  " << (coverage_ok ? "PASS" : "FAIL") << ": flood-fill reaches comparable coverage to the full independent solve" << endl;
	cout << "  " << (agreement_ok ? "PASS" : "FAIL") << ": flood-filled displacements match the independently-solved ground truth" << endl;
	if (!coverage_ok) failures++;
	if (!agreement_ok) failures++;

	//--- regression: a POI reached first by a marginal neighbor must still be solvable
	//by a better neighbor discovered later, not permanently stranded ---
	cout << endl << "=== Regression: a failed attempt from one neighbor doesn't strand the cell for a better one ===" << endl;
	{
		//1D row of 4 POIs at x = 0,1,2,3 (idx0..idx3). idx0 ("A") is a seed with high
		//ZNCC but a deliberately WRONG gradient (ux=-0.5, true is +0.5) -- realistic in
		//that a solver's own converged ZNCC doesn't guarantee its local gradient
		//estimate is equally reliable. idx3 ("B0") is a seed with lower ZNCC but the
		//CORRECT gradient. idx1 ("C") is adjacent to both A directly, and to B0's own
		//neighbor idx2 ("B") -- reachable from either side.
		//
		//Pop order is strictly ZNCC-descending, so A (higher ZNCC) is tried before B0
		//gets any chance to reach C: A's wrong gradient sends the stub a bad guess for
		//C, which fails. B0 then succeeds at B, and B then reaches C with a correct
		//guess. Whether C ends up correctly solved depends entirely on whether A's
		//earlier failed attempt was allowed to permanently claim it.
		vector<POI2D> chain;
		for (int x = 0; x < 4; x++) chain.push_back(POI2D(Point2D((float)x, 0.f)));

		chain[0].deformation.u = 0.f; //true value at x=0 (coincidentally correct)
		chain[0].deformation.ux = -0.5f; //deliberately wrong gradient
		chain[0].result.zncc = 0.99f;

		chain[3].deformation.u = 1.5f; //true value at x=3
		chain[3].deformation.ux = 0.5f; //correct gradient
		chain[3].result.zncc = 0.7f;

		LinearFieldStubSolver stub;
		ReliabilityGuided2D rg_race(4, 1);
		rg_race.setZnccThreshold(0.5f);
		rg_race.setDeltaDispTolerance(1000.f); //not the mechanism under test here
		rg_race.compute(chain, { 0, 3 }, stub);

		cout << "  C (idx1): u=" << chain[1].deformation.u << " zncc=" << chain[1].result.zncc << endl;
		bool race_ok = chain[1].result.zncc == 0.9f && std::fabs(chain[1].deformation.u - 0.5f) < 1e-4f;
		cout << "  " << (race_ok ? "PASS" : "FAIL")
			<< ": a later, better neighbor (via B0->B) still solves C after A's earlier attempt failed" << endl;
		if (!race_ok) failures++;
	}

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;
	return failures == 0 ? 0 : 1;
}
