/*
 * This file is part of OpenCorr, an open source C++ library for
 * study and development of 2D, 3D/stereo and volumetric
 * digital image correlation.
 *
 * Copyright (C) 2021-2025, Zhenyu Jiang <zhenyujiang@scut.edu.cn>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one from http://mozilla.org/MPL/2.0/.
 *
 * More information about OpenCorr can be found at https://www.opencorr.org/
 */

#include <fstream>
#include <iomanip>
#include <limits>

#include "oc_io.h"

namespace opencorr
{
	namespace
	{
		//splits a CSV data line into floats, skipping blank fields -- shared by every
		//load*() function below (loadPoint2D/3D, loadTable2D/3D, loadDeformationTable2D,
		//loadTable2DS), previously six separate copies of the same tokenizing loop, one
		//of which (loadPoint3D's own hand-rolled 3-field version) hadn't received the
		//same npos-wraparound fix already applied to loadPoint2D: a missing delimiter on
		//its second field fed a huge position2-position1 length into substr() (relying on
		//substr() to clamp it, rather than failing cleanly). Using this shared, already-
		//robust version everywhere closes that gap too, not just the duplication itself.
		std::vector<float> tokenizeCsvLine(const std::string& data_line, const std::string& delimiter)
		{
			std::vector<float> key_buffer;
			size_t position1 = 0, position2 = 0;
			std::string variable;

			do
			{
				position2 = data_line.find(delimiter, position1);
				if (position2 == std::string::npos)
				{
					position2 = data_line.length();
				}

				variable = data_line.substr(position1, position2 - position1);
				if (!variable.empty())
				{
					key_buffer.push_back(std::stof(variable));
				}

				position1 = position2 + delimiter.length();
			} while (position2 < data_line.length() && position1 < data_line.length());

			return key_buffer;
		}
	}

	IO2D::IO2D() {}

	IO2D::~IO2D() {}

	std::string IO2D::getPath() const
	{
		return file_path;
	}

	std::string IO2D::getDelimiter() const
	{
		return delimiter;
	}

	int IO2D::getWidth() const
	{
		return width;
	}

	int IO2D::getHeight()const
	{
		return height;
	}

	void IO2D::setPath(std::string file_path)
	{
		this->file_path = file_path;
	}

	void IO2D::setDelimiter(std::string delimiter)
	{
		this->delimiter = delimiter;
	}

	void IO2D::setWidth(int width)
	{
		this->width = width;
	}

	void IO2D::setHeight(int height)
	{
		this->height = height;
	}

	std::vector<Point2D> IO2D::loadPoint2D(std::string file_path)
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<Point2D> point_queue;
		int point_number = 0;

		while (getline(file_in, data_line))
		{
			point_number++;
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			if (key_buffer.size() < 2)
			{
				std::cerr << "skipping malformed row (too few fields) at line: " << point_number << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			Point2D current_point(x, y);

			point_queue.push_back(current_point);
		}
		file_in.close();

		return point_queue;
	}

	void IO2D::savePoint2D(std::vector<Point2D> point_queue, std::string file_path)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(4);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;
			file_out << std::endl;

			for (std::vector<Point2D>::iterator iter = point_queue.begin(); iter != point_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	void IO2D::loadCalibration(Calibration& calibration_cam1, Calibration& calibration_cam2, std::string file_path)
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		size_t position1, position2;
		std::string variable;

		//read intrinsics
		for (int i = 0; i < 13; i++)
		{
			getline(file_in, data_line);
			position1 = 0;
			position2 = 0;
			position2 = data_line.find(delimiter, position1);
			if (position2 == std::string::npos)
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}

			//read parameter for cam1
			position1 = position2 + delimiter.length();
			position2 = data_line.find(delimiter, position1);
			if (position2 == std::string::npos)
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}
			else
			{
				variable = data_line.substr(position1, position2 - position1);
				if (!variable.empty())
				{
					calibration_cam1.intrinsics.cam_i[i] = std::stof(variable);
				}
				else
				{
					std::cerr << "failed to read parameter at line: " << i << std::endl;
				}
			}

			//read parameter for cam2
			position1 = position2 + delimiter.length();
			position2 = data_line.length();
			variable = data_line.substr(position1, position2 - position1);
			if (!variable.empty())
			{
				calibration_cam2.intrinsics.cam_i[i] = std::stof(variable);
			}
			else
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}
		}

		//read extrinsics
		for (int i = 0; i < 6; i++)
		{
			getline(file_in, data_line);
			position1 = 0;
			position2 = 0;
			position2 = data_line.find(delimiter, position1);
			if (position2 == std::string::npos)
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}

			//read parameter for cam1
			position1 = position2 + delimiter.length();
			position2 = data_line.find(delimiter, position1);
			if (position2 == std::string::npos)
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}
			else
			{
				variable = data_line.substr(position1, position2 - position1);
				if (!variable.empty())
				{
					calibration_cam1.extrinsics.cam_e[i] = std::stof(variable);
				}
				else
				{
					std::cerr << "failed to read parameter at line: " << i << std::endl;
				}
			}

			//read parameter for cam2
			position1 = position2 + delimiter.length();
			position2 = data_line.length();
			variable = data_line.substr(position1, position2 - position1);
			if (!variable.empty())
			{
				calibration_cam2.extrinsics.cam_e[i] = std::stof(variable);
			}
			else
			{
				std::cerr << "failed to read parameter at line: " << i << std::endl;
			}
		}
		file_in.close();
	}

	std::vector<POI2D> IO2D::loadTable2D()
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<POI2D> poi_queue;

		while (getline(file_in, data_line))
		{
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			//a blank line, a row with a missing/stray delimiter, or a truncated file previously
			//indexed key_buffer[0]/[1]/etc. with no size check at all -- undefined behavior
			//(garbage values, or a crash) instead of a clean, reported skip
			if (key_buffer.size() < 4)
			{
				std::cerr << "skipping malformed row (too few fields): " << data_line << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			POI2D current_POI(x, y);

			current_POI.deformation.u = key_buffer[2];
			current_POI.deformation.v = key_buffer[3];

			int current_index = 4;
			int max_result_size = (int)(sizeof(current_POI.result.r) / sizeof(current_POI.result.r[0]));
			int strain_size = (int)(sizeof(current_POI.strain.e) / sizeof(current_POI.strain.e[0]));
			//infer the actual result-column count from the row width rather than assuming
			//the current struct size: a legacy file saved before sigma/beta were added has
			//two fewer result columns (6 vs 8), and reading it against the current 8-float
			//size would misread strain/subset_radius from the wrong offsets entirely
			int array_size = (int)key_buffer.size() - current_index - strain_size - 2;
			if (array_size > max_result_size) array_size = max_result_size;
			if (array_size < 0) array_size = 0;

			//the clamp above only protects the result.r loop itself -- it doesn't guarantee
			//enough fields remain for the strain/subset_radius reads that follow, which
			//previously ran unconditionally regardless of how short (or how mis-clamped) the
			//row actually was
			if ((int)key_buffer.size() < current_index + array_size + strain_size + 2)
			{
				std::cerr << "skipping malformed/truncated row (fewer fields than its own inferred width implies): " << data_line << std::endl;
				continue;
			}

			for (int i = 0; i < array_size; i++)
			{
				current_POI.result.r[i] = key_buffer[current_index + i];
			}
			if (array_size < max_result_size)
			{
				//legacy file predates sigma/beta -- mark them as not computed rather than
				//leaving whatever POI2D::clear() happened to zero-initialize
				current_POI.result.sigma = -1.f;
				current_POI.result.beta = 0.f;
			}

			current_index += array_size;
			array_size = (int)(sizeof(current_POI.strain.e) / sizeof(current_POI.strain.e[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.strain.e[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.subset_radius.x = key_buffer[current_index];
			current_POI.subset_radius.y = key_buffer[current_index + 1];

			poi_queue.push_back(current_POI);
		}
		file_in.close();

		return poi_queue;
	}

	void IO2D::saveTable2D(std::vector<POI2D>& poi_queue)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;

			file_out << "u" << delimiter;
			file_out << "v" << delimiter;

			file_out << "u0" << delimiter;
			file_out << "v0" << delimiter;
			file_out << "ZNCC" << delimiter;
			file_out << "iteration" << delimiter;
			file_out << "convergence" << delimiter;
			file_out << "feature" << delimiter;
			file_out << "sigma" << delimiter;
			file_out << "beta" << delimiter;

			file_out << "exx" << delimiter;
			file_out << "eyy" << delimiter;
			file_out << "exy" << delimiter;

			file_out << "subset_rx" << delimiter;
			file_out << "subset_ry" << delimiter;
			file_out << std::endl;

			for (std::vector<POI2D>::iterator iter = poi_queue.begin(); iter != poi_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;

				file_out << iter->deformation.u << delimiter;
				file_out << iter->deformation.v << delimiter;

				int array_size = (int)(sizeof(iter->result.r) / sizeof(iter->result.r[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->result.r[i] << delimiter;
				}

				array_size = (int)(sizeof(iter->strain.e) / sizeof(iter->strain.e[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->strain.e[i] << delimiter;
				}

				file_out << iter->subset_radius.x << delimiter;
				file_out << iter->subset_radius.y << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	void IO2D::saveDeformationTable2D(std::vector<POI2D>& poi_queue)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;

			file_out << "u" << delimiter;
			file_out << "ux" << delimiter;
			file_out << "uy" << delimiter;
			file_out << "uxx" << delimiter;
			file_out << "uxy" << delimiter;
			file_out << "uyy" << delimiter;

			file_out << "v" << delimiter;
			file_out << "vx" << delimiter;
			file_out << "vy" << delimiter;
			file_out << "vxx" << delimiter;
			file_out << "vxy" << delimiter;
			file_out << "vyy" << delimiter;

			file_out << "subset_rx" << delimiter;
			file_out << "subset_ry" << delimiter;
			file_out << std::endl;

			for (std::vector<POI2D>::iterator iter = poi_queue.begin(); iter != poi_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;

				int array_size = (int)(sizeof(iter->deformation.p) / sizeof(iter->deformation.p[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->deformation.p[i] << delimiter;
				}

				file_out << iter->subset_radius.x << delimiter;
				file_out << iter->subset_radius.y << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	std::vector<POI2D> IO2D::loadDeformationTable2D()
	{
		//counterpart to saveDeformationTable2D() (see oc_io.h) -- reads back x, y, the full
		//12-component deformation vector, and subset_radius, matching that function's own
		//column layout exactly
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<POI2D> poi_queue;

		while (getline(file_in, data_line))
		{
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			POI2D size_probe(0.f, 0.f);
			size_t min_fields = 2 //x, y
				+ sizeof(size_probe.deformation.p) / sizeof(size_probe.deformation.p[0])
				+ 2; //subset_radius
			if (key_buffer.size() < min_fields)
			{
				std::cerr << "skipping malformed row (too few fields): " << data_line << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			POI2D current_POI(x, y);

			int current_index = 2;
			int array_size = (int)(sizeof(current_POI.deformation.p) / sizeof(current_POI.deformation.p[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.deformation.p[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.subset_radius.x = key_buffer[current_index];
			current_POI.subset_radius.y = key_buffer[current_index + 1];

			poi_queue.push_back(current_POI);
		}
		file_in.close();

		return poi_queue;
	}

	void IO2D::saveMap2D(std::vector<POI2D>& poi_queue, OutputVariable variable)
	{
		int height = getHeight();
		int width = getWidth();
		Eigen::MatrixXf output_map = Eigen::MatrixXf::Zero(height, width);

		switch (variable)
		{
		case u:
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.u;
			}
			break;
		case v:
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.v;
			}
			break;
		case zncc: //ZNCC value
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.zncc;
			}
			break;
		case deformation_increment: //final ||delta_p||
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.convergence;
			}
			break;
		case iteration_step: //iteration steps
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.iteration;
			}
			break;
		case feature_nearby: //number of neighbor features
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.feature;
			}
			break;
		case e_xx: //strain exx
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.exx;
			}
			break;
		case e_yy: //strain eyy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.eyy;
			}
			break;
		case e_xy: //strain exy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.exy;
			}
			break;
		case u_x: //displacement gradient du/dx
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.ux;
			}
			break;
		case u_y: //displacement gradient du/dy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.uy;
			}
			break;
		case v_x: //displacement gradient dv/dx
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.vx;
			}
			break;
		case v_y: //displacement gradient dv/dy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.vy;
			}
			break;
		default:
			return;
		}

		//a failed POI's result.zncc is one of the solvers' own StatusFlag failure sentinels
		//(oc_dic.h), not a real correlation value -- if the correlation itself failed, NONE of
		//that POI's fields are trustworthy, whichever `variable` was requested. Overwriting
		//with NaN (rather than leaving the raw sentinel, e.g. -8) matters specifically for this
		//map/raster export: a viewer (ParaView, or any tool coloring this as a heatmap) would
		//otherwise render a failure as if it were a real, very-poor-but-genuine correlation,
		//skewing any color scale or statistic computed over the map. NaN is the raster
		//convention such tools already treat as "no data" -- 0 isn't safe here since it's
		//already this map's own "no POI was ever placed at this pixel" default.
		for (size_t i = 0; i < poi_queue.size(); i++)
		{
			if (poi_queue[i].result.zncc < 0.f)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = std::numeric_limits<float>::quiet_NaN();
			}
		}

		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);
		if (file_out.is_open())
		{
			for (int r = 0; r < height; r++)
			{
				for (int c = 0; c < width; c++)
				{
					file_out << output_map(r, c) << delimiter;
				}
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	std::vector<POI2DS> IO2D::loadTable2DS()
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<POI2DS> poi_queue;

		while (getline(file_in, data_line))
		{
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			POI2DS size_probe(0.f, 0.f);
			size_t min_fields = 2
				+ sizeof(size_probe.deformation.p) / sizeof(size_probe.deformation.p[0])
				+ sizeof(size_probe.result.r) / sizeof(size_probe.result.r[0])
				+ 6 //ref_coor + tar_coor
				+ sizeof(size_probe.strain.e) / sizeof(size_probe.strain.e[0])
				+ 2; //subset_radius
			if (key_buffer.size() < min_fields)
			{
				std::cerr << "skipping malformed row (too few fields): " << data_line << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			POI2DS current_POI(x, y);

			int current_index = 2;
			int array_size = (int)(sizeof(current_POI.deformation.p) / sizeof(current_POI.deformation.p[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.deformation.p[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			array_size = (int)(sizeof(current_POI.result.r) / sizeof(current_POI.result.r[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.result.r[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.ref_coor.x = key_buffer[current_index];
			current_POI.ref_coor.y = key_buffer[current_index + 1];
			current_POI.ref_coor.z = key_buffer[current_index + 2];

			current_POI.tar_coor.x = key_buffer[current_index + 3];
			current_POI.tar_coor.y = key_buffer[current_index + 4];
			current_POI.tar_coor.z = key_buffer[current_index + 5];

			current_index += 6;
			array_size = (int)(sizeof(current_POI.strain.e) / sizeof(current_POI.strain.e[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.strain.e[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.subset_radius.x = key_buffer[current_index];
			current_POI.subset_radius.y = key_buffer[current_index + 1];

			poi_queue.push_back(current_POI);
		}
		file_in.close();

		return poi_queue;
	}

	void IO2D::saveTable2DS(std::vector<POI2DS>& poi_queue)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;

			file_out << "u" << delimiter;
			file_out << "v" << delimiter;
			file_out << "w" << delimiter;

			file_out << "r1r2 ZNCC" << delimiter;
			file_out << "r1t1 ZNCC" << delimiter;
			file_out << "r1t2 ZNCC" << delimiter;

			file_out << "r2_x" << delimiter;
			file_out << "r2_y" << delimiter;
			file_out << "t1_x" << delimiter;
			file_out << "t1_y" << delimiter;
			file_out << "t2_x" << delimiter;
			file_out << "t2_y" << delimiter;

			file_out << "ref_x" << delimiter;
			file_out << "ref_y" << delimiter;
			file_out << "ref_z" << delimiter;
			file_out << "tar_x" << delimiter;
			file_out << "tar_y" << delimiter;
			file_out << "tar_z" << delimiter;

			file_out << "exx" << delimiter;
			file_out << "eyy" << delimiter;
			file_out << "ezz" << delimiter;
			file_out << "exy" << delimiter;
			file_out << "eyz" << delimiter;
			file_out << "ezx" << delimiter;

			file_out << "subset_rx" << delimiter;
			file_out << "subset_ry" << delimiter;
			file_out << std::endl;

			for (std::vector<POI2DS>::iterator iter = poi_queue.begin(); iter != poi_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;

				int array_size = (int)(sizeof(iter->deformation.p) / sizeof(iter->deformation.p[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->deformation.p[i] << delimiter;
				}

				array_size = (int)(sizeof(iter->result.r) / sizeof(iter->result.r[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->result.r[i] << delimiter;
				}

				file_out << iter->ref_coor.x << delimiter;
				file_out << iter->ref_coor.y << delimiter;
				file_out << iter->ref_coor.z << delimiter;

				file_out << iter->tar_coor.x << delimiter;
				file_out << iter->tar_coor.y << delimiter;
				file_out << iter->tar_coor.z << delimiter;

				array_size = (int)(sizeof(iter->strain.e) / sizeof(iter->strain.e[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->strain.e[i] << delimiter;
				}

				file_out << iter->subset_radius.x << delimiter;
				file_out << iter->subset_radius.y << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	void IO2D::saveMap2DS(std::vector<POI2DS>& poi_queue, OutputVariable  variable)
	{
		int height = getHeight();
		int width = getWidth();
		Eigen::MatrixXf output_map = Eigen::MatrixXf::Zero(height, width);

		switch (variable)
		{
		case u:
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.u;
			}
			break;
		case v:
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.v;
			}
			break;
		case w:
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].deformation.w;
			}
			break;
		case zncc_r1r2: //ZNCC value in matching between the two reference images
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.r1r2_zncc;
			}
			break;
		case zncc: //ZNCC value in matching between the reference image and the target image from same view
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.r1t1_zncc;
			}
			break;
		case zncc_r1t2: //ZNCC value in matching between the reference image and the target image from different views
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].result.r1t2_zncc;
			}
			break;
		case e_xx: //strain exx
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.exx;
			}
			break;
		case e_yy: //strain eyy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.eyy;
			}
			break;
		case e_zz: //strain ezz
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.ezz;
			}
			break;
		case e_xy: //strain exy
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.exy;
			}
			break;
		case e_yz: //strain eyz
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.eyz;
			}
			break;
		case e_zx: //strain ezx
			for (int i = 0; i < poi_queue.size(); i++)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = poi_queue[i].strain.ezx;
			}
			break;
		default:
			return;
		}

		//see saveMap2D's own comment on why this matters -- a POI2DS carries three
		//independent ZNCC values (one per matching stage: r1-r2 stereo, r1-t1 temporal,
		//r1-t2 stereo); if ANY of them failed, none of this POI's fields (whichever
		//`variable` was requested) are trustworthy, so the whole cell becomes NaN
		for (size_t i = 0; i < poi_queue.size(); i++)
		{
			if (poi_queue[i].result.r1r2_zncc < 0.f || poi_queue[i].result.r1t1_zncc < 0.f || poi_queue[i].result.r1t2_zncc < 0.f)
			{
				output_map((int)poi_queue[i].y, (int)poi_queue[i].x) = std::numeric_limits<float>::quiet_NaN();
			}
		}

		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			for (int r = 0; r < height; r++)
			{
				for (int c = 0; c < width; c++)
				{
					file_out << output_map(r, c) << delimiter;
				}
				file_out << std::endl;
			}
		}
		file_out.close();
	}


	IO3D::IO3D() {}

	IO3D::~IO3D() {}

	std::string IO3D::getPath() const
	{
		return file_path;
	}

	std::string IO3D::getDelimiter() const
	{
		return delimiter;
	}

	void IO3D::setPath(std::string file_path)
	{
		this->file_path = file_path;
	}

	void IO3D::setDelimiter(std::string delimiter)
	{
		this->delimiter = delimiter;
	}

	int IO3D::getDimX()
	{
		return dim_x;
	}

	int IO3D::getDimY()
	{
		return dim_y;
	}

	int IO3D::getDimZ()
	{
		return dim_z;
	}

	void IO3D::setDimX(int dim_x)
	{
		this->dim_x = dim_x;
	}

	void IO3D::setDimY(int dim_y)
	{
		this->dim_y = dim_y;
	}

	void IO3D::setDimZ(int dim_z)
	{
		this->dim_z = dim_z;
	}

	std::vector<Point3D> IO3D::loadPoint3D(std::string file_path)
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<Point3D> point_queue;
		int point_number = 0;

		while (getline(file_in, data_line))
		{
			point_number++;
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			if (key_buffer.size() < 3)
			{
				std::cerr << "skipping malformed row (too few fields) at line: " << point_number << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			float z = key_buffer[2];
			Point3D current_point(x, y, z);

			point_queue.push_back(current_point);
		}
		file_in.close();

		return point_queue;
	}

	void IO3D::savePoint3D(std::vector<Point3D> poi_queue, std::string file_path)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(4);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;
			file_out << "z" << delimiter;
			file_out << std::endl;

			for (std::vector<Point3D>::iterator iter = poi_queue.begin(); iter != poi_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;
				file_out << iter->z << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	std::vector<POI3D> IO3D::loadTable3D()
	{
		std::ifstream file_in(file_path);
		if (!file_in.is_open())
		{
			std::cerr << "failed to read file " << file_path << std::endl;
		}

		std::string data_line;
		getline(file_in, data_line);
		std::vector<POI3D> poi_queue;

		while (getline(file_in, data_line))
		{
			std::vector<float> key_buffer = tokenizeCsvLine(data_line, delimiter);

			POI3D size_probe(0.f, 0.f, 0.f);
			size_t min_fields = 3 //x, y, z
				+ 3 //u, v, w
				+ sizeof(size_probe.result.r) / sizeof(size_probe.result.r[0])
				+ 9 //ux, uy, uz, vx, vy, vz, wx, wy, wz
				+ sizeof(size_probe.strain.e) / sizeof(size_probe.strain.e[0])
				+ 3; //subset_radius
			if (key_buffer.size() < min_fields)
			{
				std::cerr << "skipping malformed row (too few fields): " << data_line << std::endl;
				continue;
			}

			float x = key_buffer[0];
			float y = key_buffer[1];
			float z = key_buffer[2];
			POI3D current_POI(x, y, z);

			current_POI.deformation.u = key_buffer[3];
			current_POI.deformation.v = key_buffer[4];
			current_POI.deformation.w = key_buffer[5];

			int current_index = 6;
			int array_size = (int)(sizeof(current_POI.result.r) / sizeof(current_POI.result.r[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.result.r[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.deformation.ux = key_buffer[current_index];
			current_POI.deformation.uy = key_buffer[current_index + 1];
			current_POI.deformation.uz = key_buffer[current_index + 2];
			current_POI.deformation.vx = key_buffer[current_index + 3];
			current_POI.deformation.vy = key_buffer[current_index + 4];
			current_POI.deformation.vz = key_buffer[current_index + 5];
			current_POI.deformation.wx = key_buffer[current_index + 6];
			current_POI.deformation.wy = key_buffer[current_index + 7];
			current_POI.deformation.wz = key_buffer[current_index + 8];
			array_size = 9;

			current_index += array_size;
			array_size = (int)(sizeof(current_POI.strain.e) / sizeof(current_POI.strain.e[0]));
			for (int i = 0; i < array_size; i++)
			{
				current_POI.strain.e[i] = key_buffer[current_index + i];
			}

			current_index += array_size;
			current_POI.subset_radius.x = key_buffer[current_index];
			current_POI.subset_radius.y = key_buffer[current_index + 1];
			current_POI.subset_radius.z = key_buffer[current_index + 2];

			poi_queue.push_back(current_POI);
		}
		file_in.close();

		return poi_queue;
	}

	void IO3D::saveTable3D(std::vector<POI3D>& poi_queue)
	{
		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			file_out << "x" << delimiter;
			file_out << "y" << delimiter;
			file_out << "z" << delimiter;

			file_out << "u" << delimiter;
			file_out << "v" << delimiter;
			file_out << "w" << delimiter;

			file_out << "u0" << delimiter;
			file_out << "v0" << delimiter;
			file_out << "w0" << delimiter;
			file_out << "ZNCC" << delimiter;
			file_out << "iteration" << delimiter;
			file_out << "convergence" << delimiter;
			file_out << "feature" << delimiter;

			file_out << "ux" << delimiter;
			file_out << "uy" << delimiter;
			file_out << "uz" << delimiter;
			file_out << "vx" << delimiter;
			file_out << "vy" << delimiter;
			file_out << "vz" << delimiter;
			file_out << "wx" << delimiter;
			file_out << "wy" << delimiter;
			file_out << "wz" << delimiter;

			file_out << "exx" << delimiter;
			file_out << "eyy" << delimiter;
			file_out << "ezz" << delimiter;
			file_out << "exy" << delimiter;
			file_out << "eyz" << delimiter;
			file_out << "ezx" << delimiter;

			file_out << "subset_rx" << delimiter;
			file_out << "subset_ry" << delimiter;
			file_out << "subset_rz" << delimiter;
			file_out << std::endl;

			for (std::vector<POI3D>::iterator iter = poi_queue.begin(); iter != poi_queue.end(); iter++)
			{
				file_out << iter->x << delimiter;
				file_out << iter->y << delimiter;
				file_out << iter->z << delimiter;

				file_out << iter->deformation.u << delimiter;
				file_out << iter->deformation.v << delimiter;
				file_out << iter->deformation.w << delimiter;

				int array_size = (int)(sizeof(iter->result.r) / sizeof(iter->result.r[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->result.r[i] << delimiter;
				}

				file_out << iter->deformation.ux << delimiter;
				file_out << iter->deformation.uy << delimiter;
				file_out << iter->deformation.uz << delimiter;
				file_out << iter->deformation.vx << delimiter;
				file_out << iter->deformation.vy << delimiter;
				file_out << iter->deformation.vz << delimiter;
				file_out << iter->deformation.wx << delimiter;
				file_out << iter->deformation.wy << delimiter;
				file_out << iter->deformation.wz << delimiter;

				array_size = (int)(sizeof(iter->strain.e) / sizeof(iter->strain.e[0]));
				for (int i = 0; i < array_size; i++)
				{
					file_out << iter->strain.e[i] << delimiter;
				}

				file_out << iter->subset_radius.x << delimiter;
				file_out << iter->subset_radius.y << delimiter;
				file_out << iter->subset_radius.z << delimiter;
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	void IO3D::saveMap3D(std::vector<POI3D>& poi_queue, OutputVariable variable)
	{
		int queue_length = (int)poi_queue.size();
		float*** output_map = new3D(getDimZ(), getDimY(), getDimX());

		switch (variable)
		{
		case u:
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.u;
			}
			break;
		case v:
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.v;
			}
			break;
		case w:
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.w;
			}
			break;
		case zncc: //ZNCC value
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].result.zncc;
			}
			break;
		case iteration_step: //iteration step
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].result.iteration;
			}
			break;
		case deformation_increment: //deformation increment
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].result.convergence;
			}
			break;
		case feature_nearby: //features nearby
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].result.feature;
			}
			break;
		case e_xx: //strain exx
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.exx;
			}
			break;
		case e_yy: //strain eyy
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.eyy;
			}
			break;
		case e_zz: //strain ezz
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.ezz;
			}
			break;
		case e_xy: //strain exy
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.exy;
			}
			break;
		case e_yz: //strain eyz
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.eyz;
			}
			break;
		case e_zx: //strain ezx
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].strain.ezx;
			}
			break;
		case u_x: //displacement gradient du/dx
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.ux;
			}
			break;
		case u_y: //displacement gradient du/dy
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.uy;
			}
			break;
		case u_z: //displacement gradient du/dz
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.uz;
			}
			break;
		case v_x: //displacement gradient dv/dx
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.vx;
			}
			break;
		case v_y: //displacement gradient dv/dy
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.vy;
			}
			break;
		case v_z: //displacement gradient dv/dz
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.vz;
			}
			break;
		case w_x: //displacement gradient dw/dx
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.wx;
			}
			break;
		case w_y: //displacement gradient dw/dy
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.wy;
			}
			break;
		case w_z: //displacement gradient dw/dz
			for (int i = 0; i < queue_length; i++)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = poi_queue[i].deformation.wz;
			}
			break;
		default:
			return;
		}

		//see IO2D::saveMap2D's own comment on why this matters
		for (int i = 0; i < queue_length; i++)
		{
			if (poi_queue[i].result.zncc < 0.f)
			{
				output_map[(int)poi_queue[i].z][(int)poi_queue[i].y][(int)poi_queue[i].x] = std::numeric_limits<float>::quiet_NaN();
			}
		}

		std::ofstream file_out(file_path);
		file_out.setf(std::ios::fixed);
		file_out << std::setprecision(8);

		if (file_out.is_open())
		{
			for (int i = 0; i < getDimZ(); i++)
			{
				for (int j = 0; j < getDimY(); j++)
				{
					for (int k = 0; k < getDimX(); k++)
					{
						file_out << output_map[i][j][k] << delimiter;
					}
					file_out << std::endl;
				}
				file_out << std::endl;
			}
		}
		file_out.close();
	}

	void IO3D::saveMatrixBin(std::vector<POI3D>& poi_queue)
	{
		std::ofstream file_out;
		file_out.open(file_path, std::ios::out | std::ios::binary);

		if (!file_out.is_open())
		{
			std::cerr << "failed to open file " << file_path << std::endl;
		}

		//head information, including the number of POIs and the three dimensions of image
		int queue_length = (int)poi_queue.size();
		int head_info[4];
		head_info[0] = queue_length;
		head_info[1] = dim_x;
		head_info[2] = dim_y;
		head_info[3] = dim_z;

		//create a 1D array and fill it with output data of POIs
		int result_length = 8;
		int data_length = result_length * queue_length;
		float* data_array = new float[data_length];

		for (int i = 0; i < queue_length; i++)
		{
			data_array[i * result_length] = poi_queue[i].x;
			data_array[i * result_length + 1] = poi_queue[i].y;
			data_array[i * result_length + 2] = poi_queue[i].z;
			data_array[i * result_length + 3] = poi_queue[i].deformation.u;
			data_array[i * result_length + 4] = poi_queue[i].deformation.v;
			data_array[i * result_length + 5] = poi_queue[i].deformation.w;
			data_array[i * result_length + 6] = poi_queue[i].result.zncc;
			data_array[i * result_length + 7] = poi_queue[i].result.convergence;
		}

		//write head information
		file_out.write((char*)head_info, sizeof(head_info[0]) * 4);

		//write data
		file_out.write((char*)data_array, sizeof(data_array[0]) * data_length);

		file_out.close();

		delete[] data_array;
	}

	std::vector<POI3D> IO3D::loadMatrixBin()
	{
		std::vector<POI3D> poi_queue;

		//opened in binary mode to match saveMatrixBin's own std::ios::binary write -- opening
		//in text mode (as this previously did) performs CRLF/LF translation on the raw
		//float/int bytes on Windows, corrupting any byte that happens to match a newline code,
		//silently returning wrong coordinates/displacements
		std::ifstream file_in(file_path, std::ios::in | std::ios::binary);
		if (!file_in.is_open())
		{
			std::cerr << "failed to open file " << file_path << std::endl;
			return poi_queue; //previously fell through to read() on a closed stream instead
		}

		//read head information, validating the read actually succeeded and queue_length is
		//sane before trusting it to size an allocation -- a missing/truncated file previously
		//read head_info from an unread/garbage buffer, and a resulting negative or absurdly
		//large queue_length fed straight into `new float[data_length]` (std::bad_array_new_length
		//or an out-of-memory crash instead of a clean, empty result)
		int head_info[4];
		file_in.read((char*)head_info, sizeof(head_info[0]) * 4);
		if (!file_in.good() || head_info[0] < 0)
		{
			std::cerr << "failed to read a valid header from file " << file_path << std::endl;
			return poi_queue;
		}
		int queue_length = head_info[0];
		setDimX(head_info[1]);
		setDimY(head_info[2]);
		setDimZ(head_info[3]);

		int result_length = 8;
		int data_length = result_length * queue_length;
		std::vector<float> data_array(data_length);
		file_in.read((char*)data_array.data(), sizeof(float) * data_length);
		if ((size_t)file_in.gcount() != sizeof(float) * data_length)
		{
			std::cerr << "file " << file_path << " is truncated relative to its own header" << std::endl;
			return poi_queue;
		}
		file_in.close();

		POI3D empty_poi(0, 0, 0);
		poi_queue.resize(queue_length, empty_poi);

		for (int i = 0; i < queue_length; i++)
		{
			poi_queue[i].x = data_array[i * result_length];
			poi_queue[i].y = data_array[i * result_length + 1];
			poi_queue[i].z = data_array[i * result_length + 2];
			poi_queue[i].deformation.u = data_array[i * result_length + 3];
			poi_queue[i].deformation.v = data_array[i * result_length + 4];
			poi_queue[i].deformation.w = data_array[i * result_length + 5];
			poi_queue[i].result.zncc = data_array[i * result_length + 6];
			poi_queue[i].result.convergence = data_array[i * result_length + 7];
		}

		return poi_queue;
	}

}//namespace opencorr