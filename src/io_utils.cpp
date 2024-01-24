#include "structs.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include "sal-test-config.h"

void read_run_config(const std::string file_name, RunConfig &run_information) {
  // reads run information of file_name
  std::ifstream config_file(file_name);
  if (config_file.fail()) {
      std::cout << "namelist file at " << file_name << std::endl;
      throw std::runtime_error("namelist not found, try running from inside /bin");
  }
  std::string line, word1, word2;

  while (true) {
    getline(config_file, line);
    std::stringstream str1(line);
    getline(str1, word1, '=');
    getline(str1, word2);
    if (word1 == "use_fast") {
      if (stoi(word2) == 1) {
        run_information.use_fast = true;
      }
    } else if (word1 == "out_path") {
      run_information.out_path = word2;
    } else if (word1 == "write_output") {
      if (stoi(word2) == 1) {
        run_information.write_output = true;
      }
    } else if (word1 == "write_precision") {
      run_information.write_precision = stoi(word2);
    } else if (word1 == "radius") {
      run_information.radius = stod(word2);
    } else if (word1 == "point_count") {
      run_information.point_count = stoi(word2);
    } else if (word1 == "fast_sum_cluster_thresh") {
      run_information.fast_sum_cluster_thresh = stoi(word2);
    } else if (word1 == "theta") {
      run_information.fast_sum_theta = stod(word2);
    } else if (word1 == "fast_sum_rotate") {
      if (stoi(word2) == 0) {
        run_information.fast_sum_rotate = false;
      }
    } else if (word1 == "time") {
      run_information.time = stoi(word2);
    } else if (word1 == "degree") {
      run_information.interp_degree = stoi(word2);
    } else if (word1 == "fast_rot_alph") {
      run_information.fast_sum_rotate_alph = stod(word2);
    } else if (word1 == "fast_rot_beta") {
      run_information.fast_sum_rotate_beta = stod(word2);
    } else if (word1 == "fast_rot_gamm") {
      run_information.fast_sum_rotate_gamm = stod(word2);
    } else if (word1 == "sph_harm_comps") {
      run_information.sph_harm_comps = stoi(word2);
    } else if (word1 == "use_ces") {
      if (stoi(word2) == 1) {
        run_information.use_cesaro = true;
      }
    } else {
      run_information.interp_point_count = (run_information.interp_degree + 1) * (run_information.interp_degree + 2) / 2;
      return;
    }
  }
}

std::string create_config(const RunConfig &run_information) {
  std::stringstream ss1, ss2, ss3;
  int precision;
  std::string output_filename = std::to_string(run_information.time) + "_" +
      std::to_string(run_information.point_count) + "_" + std::to_string(run_information.sph_harm_comps) + "_";
  if (run_information.use_fast) {
    output_filename +=
        "fast_" + std::to_string(run_information.fast_sum_cluster_thresh) + "_" +
        std::to_string(run_information.fast_sum_theta).substr(0, 3) + "_" + std::to_string(run_information.interp_degree);
  } else
    output_filename += "direct";
  return output_filename;
}

void read_data(const RunConfig &run_information, std::vector<double>& latvals, std::vector<double>& lonvals, std::vector<double>& sshs,
               std::vector<double>& areas) {
  // reads in lats, lons, and ssh values
  std::ifstream fin1, fin2, fin3, fin4;
  fin1.open(DATA_DIR + std::string("mpas_grid_lats.csv"));
  fin2.open(DATA_DIR + std::string("mpas_grid_lons.csv"));
  fin3.open(DATA_DIR + std::string("mpas_grid_ssh") + std::to_string(run_information.time) + ".csv");
  fin4.open(DATA_DIR + std::string("mpas_grid_area.csv"));
  std::string line1, line2, line3, line4;
  double scaling_factor = 4 * M_PI * pow(6371229, 2);
  for (int i = 0; i < run_information.point_count; i++) {
    getline(fin1, line1, ',');
    getline(fin2, line2, ',');
    getline(fin3, line3, ',');
    getline(fin4, line4, ',');
    latvals[i] = stod(line1);
    lonvals[i] = stod(line2);
    sshs[i] = stod(line3);
    areas[i] = stod(line4) / scaling_factor;
  }
}

void load_llns(const RunConfig &run_information, std::vector<double>& llns) {
  // read the LLNs
  std::ifstream fin1, fin2;
  fin1.open(DATA_DIR + std::string("lln_hs.csv"));
  fin2.open(DATA_DIR + std::string("lln_ks.csv"));
  std::string line1, line2;
  for (int i = 0; i < run_information.sph_harm_comps; i++) {
    getline(fin1, line1, ',');
    getline(fin2, line2, ',');
    llns[i] = 1 + stod(line2) - stod(line1);
  }
}

void write_state(const RunConfig &run_information, const std::vector<double> &sal, const std::string path) {
  // write sal potential
  std::ofstream write_out(path, std::ofstream::out | std::ofstream::trunc);
  for (int i = 0; i < run_information.point_count;
       i++) { // write out initial state
    write_out << std::setprecision(run_information.write_precision) << sal[i] << "\n";
  }
  write_out.close();
}
