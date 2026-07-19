/*
 Local build-verification test for Uncertainty2D (sigma/beta), issue #2.
 Runs the same FFTCC + ICGN1 pipeline as the existing 2D DIC example, then
 computes sigma/beta on top of the converged POIs.
*/

#include <fstream>
#include <omp.h>

#include "opencorr.h"

using namespace opencorr;
using namespace std;

int main()
{
	int failures = 0;

	string ref_image_path = "examples/2d_dic/oht_cfrp_0.bmp";
	string tar_image_path = "examples/2d_dic/oht_cfrp_4.bmp";
	Image2D ref_img(ref_image_path);
	Image2D tar_img(tar_image_path);

	int cpu_thread_number = omp_get_num_procs() - 1;
	omp_set_num_threads(cpu_thread_number);

	int subset_radius_x = 16;
	int subset_radius_y = 16;
	int max_iteration = 10;
	float max_deformation_norm = 0.001f;

	Point2D upper_left_point(30, 30);
	vector<POI2D> poi_queue;
	int poi_number_x = 100;
	int poi_number_y = 300;
	int grid_space = 2;

	for (int i = 0; i < poi_number_y; i++)
	{
		for (int j = 0; j < poi_number_x; j++)
		{
			Point2D offset(j * grid_space, i * grid_space);
			POI2D current_poi(upper_left_point + offset);
			poi_queue.push_back(current_poi);
		}
	}

	cout << poi_queue.size() << " POIs, " << cpu_thread_number << " CPU threads." << endl;

	FFTCC2D* fftcc = new FFTCC2D(subset_radius_x, subset_radius_y, cpu_thread_number);
	fftcc->setImages(ref_img, tar_img);
	fftcc->compute(poi_queue);

	ICGN2D1* icgn1 = new ICGN2D1(subset_radius_x, subset_radius_y, max_deformation_norm, max_iteration, cpu_thread_number);
	icgn1->setImages(ref_img, tar_img);
	icgn1->prepare();
	icgn1->compute(poi_queue);

	int converged = 0;
	for (auto& poi : poi_queue) if (poi.result.zncc > 0.9f) converged++;
	cout << "ICGN converged (ZNCC > 0.9) on " << converged << " / " << poi_queue.size() << " POIs." << endl;

	double timer_tic = omp_get_wtime();

	Uncertainty2D* uq = new Uncertainty2D(subset_radius_x, subset_radius_y, cpu_thread_number);
	uq->setImages(ref_img, tar_img);
	uq->prepare();
	uq->compute(poi_queue);

	double consumed_time = omp_get_wtime() - timer_tic;
	cout << "Uncertainty quantification (sigma/beta) takes " << consumed_time << " sec." << endl;

	//summary stats, and a handful of sample rows
	double sum_sigma = 0.0, sum_beta = 0.0;
	int valid_sigma = 0, valid_beta = 0;
	float min_sigma = 1e9f, max_sigma = -1e9f, min_beta = 1e9f, max_beta = -1e9f;
	for (auto& poi : poi_queue)
	{
		if (poi.result.zncc <= 0.9f) continue;
		if (poi.result.sigma >= 0.f)
		{
			sum_sigma += poi.result.sigma;
			valid_sigma++;
			min_sigma = std::min(min_sigma, poi.result.sigma);
			max_sigma = std::max(max_sigma, poi.result.sigma);
		}
		//beta == -1.f: DICe-style "undefined" sentinel (cost too flat along some axis),
		//distinct from a real (always non-negative) conditioning value
		if (poi.result.beta >= 0.f)
		{
			sum_beta += poi.result.beta;
			valid_beta++;
			min_beta = std::min(min_beta, poi.result.beta);
			max_beta = std::max(max_beta, poi.result.beta);
		}
	}

	cout << "sigma: mean=" << (valid_sigma ? sum_sigma / valid_sigma : 0.0)
		<< " min=" << min_sigma << " max=" << max_sigma
		<< " (" << valid_sigma << " valid of " << converged << " converged)" << endl;
	cout << "beta:  mean=" << (valid_beta ? sum_beta / valid_beta : 0.0)
		<< " min=" << min_beta << " max=" << max_beta
		<< " (" << valid_beta << " valid of " << converged << " converged)" << endl;

	//dump a small sample table
	string file_path = tar_image_path.substr(0, tar_image_path.find_last_of(".")) + "_uncertainty_sample.csv";
	ofstream csv_out(file_path);
	csv_out << "x,y,zncc,sigma,beta" << endl;
	int printed = 0;
	for (auto& poi : poi_queue)
	{
		if (poi.result.zncc <= 0.9f) continue;
		csv_out << poi.x << "," << poi.y << "," << poi.result.zncc << "," << poi.result.sigma << "," << poi.result.beta << endl;
		if (printed < 10)
		{
			cout << "x=" << poi.x << " y=" << poi.y << " zncc=" << poi.result.zncc
				<< " sigma=" << poi.result.sigma << " beta=" << poi.result.beta << endl;
		}
		printed++;
	}
	csv_out.close();
	cout << "Wrote " << printed << " rows to " << file_path << endl;

	//--- regression: beta's "not computed"/"degenerate subset" sentinel must be -1.f, not
	//0.f -- beta is a reciprocal-slope conditioning score where smaller is better, so 0.f
	//would misleadingly read as "extremely well-conditioned" instead of "never evaluated" ---
	cout << endl << "=== beta sentinel must be -1.f, not 0.f, for POIs never evaluated ===" << endl;
	{
		//case 1: precondition failure (zncc < 0) -- never reaches the conditioning probe at all
		POI2D poi_failed_precondition(Point2D(60.f, 60.f));
		poi_failed_precondition.result.zncc = -1.f;
		Uncertainty2D uq_precondition(subset_radius_x, subset_radius_y, 1);
		uq_precondition.setImages(ref_img, tar_img);
		uq_precondition.prepare();
		uq_precondition.compute(&poi_failed_precondition);
		bool precondition_ok = poi_failed_precondition.result.beta == -1.f;
		cout << "  zncc<0 precondition: beta=" << poi_failed_precondition.result.beta
			<< " " << (precondition_ok ? "PASS" : "FAIL") << endl;
		if (!precondition_ok) failures++;

		//case 2: out-of-bounds subset (too close to the image edge)
		POI2D poi_out_of_bounds(Point2D(2.f, 2.f)); //subset_radius=16, so this subset runs off the image
		poi_out_of_bounds.result.zncc = 0.9f;
		Uncertainty2D uq_bounds(subset_radius_x, subset_radius_y, 1);
		uq_bounds.setImages(ref_img, tar_img);
		uq_bounds.prepare();
		uq_bounds.compute(&poi_out_of_bounds);
		bool bounds_ok = poi_out_of_bounds.result.beta == -1.f;
		cout << "  out-of-bounds subset: beta=" << poi_out_of_bounds.result.beta
			<< " " << (bounds_ok ? "PASS" : "FAIL") << endl;
		if (!bounds_ok) failures++;

		//case 3: degenerate (perfectly uniform-intensity) subset -- ref_mean_norm <= 0
		Image2D flat_ref(64, 64), flat_tar(64, 64);
		flat_ref.eg_mat.setConstant(128.f);
		flat_tar.eg_mat.setConstant(128.f);
		POI2D poi_flat(Point2D(32.f, 32.f));
		poi_flat.result.zncc = 0.9f;
		Uncertainty2D uq_flat(subset_radius_x, subset_radius_y, 1);
		uq_flat.setImages(flat_ref, flat_tar);
		uq_flat.prepare();
		uq_flat.compute(&poi_flat);
		bool flat_ok = poi_flat.result.beta == -1.f;
		cout << "  degenerate uniform-intensity subset: beta=" << poi_flat.result.beta
			<< " " << (flat_ok ? "PASS" : "FAIL") << endl;
		if (!flat_ok) failures++;

		//case 4: znssd()'s own out-of-bounds guard -- a POI whose REFERENCE position passes
		//compute()'s coarse bounds precondition, but whose already-converged deformation.u
		//warps the TARGET subset entirely outside the image. Without znssd()'s own guard,
		//BicubicBspline's -1.f out-of-bounds sentinel would be sampled as if it were a real
		//pixel intensity and silently blended into gamma_0/gamma_p/gamma_m.
		POI2D poi_warped_oob(Point2D(32.f, 32.f));
		poi_warped_oob.result.zncc = 0.9f;
		poi_warped_oob.deformation.u = 1000.f; //warps the target subset far outside any real image
		Uncertainty2D uq_warped_oob(subset_radius_x, subset_radius_y, 1);
		uq_warped_oob.setImages(ref_img, tar_img);
		uq_warped_oob.prepare();
		uq_warped_oob.compute(&poi_warped_oob);
		bool warped_oob_ok = poi_warped_oob.result.beta == -1.f;
		cout << "  converged deformation warps target subset out of bounds: beta=" << poi_warped_oob.result.beta
			<< " " << (warped_oob_ok ? "PASS" : "FAIL") << endl;
		if (!warped_oob_ok) failures++;
	}

	cout << endl << (failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << " (" << failures << " failure(s))" << endl;

	delete fftcc;
	delete icgn1;
	delete uq;

	return failures == 0 ? 0 : 1;
}
