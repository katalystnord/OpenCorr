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

	delete fftcc;
	delete icgn1;
	delete uq;

	return 0;
}
