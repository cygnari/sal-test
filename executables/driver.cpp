#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "sal-test-config.h"
#include "structs.hpp"
#include "io_utils.hpp"
#include "green_func.hpp"
#include "general_utils.hpp"
#include "icos_funcs.hpp"
#include "fast_sum_funcs.hpp"
#include "mpi_utils.hpp"

int main(int argc, char **argv) {
  // initialize everything
  MPI_Init(&argc, &argv);
  int P, ID;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &ID);
  MPI_Win win_sal;

  RunConfig run_information;
  const std::string namelist_file = std::string(NAMELIST_DIR) + std::string("namelist.txt");
  read_run_config(namelist_file, run_information);

  std::cout << run_information.point_count << std::endl;

  run_information.mpi_P = P;
  run_information.mpi_ID = ID;

  std::vector<double> lats (run_information.point_count, 0); // latitudes of the points, -pi/2 to pi/2
  std::vector<double> lons (run_information.point_count, 0); // longitudes, 0 to 2pi
  std::vector<double> xval (run_information.point_count, 0); // x vals in R^3
  std::vector<double> yval (run_information.point_count, 0); // y vals in R^3
  std::vector<double> zval (run_information.point_count, 0); // z vals in R^3
  std::vector<double> sshs (run_information.point_count, 0); // sea surface height, meters
  std::vector<double> sals (run_information.point_count, 0); // SAL potential
  std::vector<double> area (run_information.point_count, 0); // area of each grid cell

  MPI_Win_create(&sals[0], run_information.point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sal);

  std::vector<double> llns (run_information.sph_harm_comps, 0); // load love numbers

  std::vector<IcosPanel> fast_sum_icos_panels; // faces of the fast sum icosahedron, tree
  std::vector<InteractPair> fast_sum_interactions; // result of dual tree traversal

  std::string output_folder = create_config(run_information);

  if (ID == 0) {
    std::string filename = NAMELIST_DIR +
        std::string("initialize.py ") + run_information.out_path + "/" + output_folder;
    std::string command = "python ";
    command += filename;
    system(command.c_str());
  }

  // read in data

  read_data(run_information, lats, lons, sshs, area);
  load_llns(run_information, llns);

  std::tuple<double, double, double> xyz;
  std::chrono::steady_clock::time_point begin, end;

  if (run_information.use_fast) {
    for (int i = 0; i < run_information.point_count; i++) {
      xyz = latlon_to_xyz(lats[i], lons[i], run_information.radius);
      xval[i] = std::get<0>(xyz);
      yval[i] = std::get<1>(xyz);
      zval[i] = std::get<2>(xyz);
    }
    if (ID == 0) {
      begin = std::chrono::steady_clock::now();
    }

    initialize_icosahedron(run_information, fast_sum_icos_panels, xval, yval, zval);

    if (ID == 0) {
      end = std::chrono::steady_clock::now();
      std::cout << "icos init time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()
                << " microseconds" << std::endl;
    }

    if (ID == 0) {
      begin = std::chrono::steady_clock::now();
    }

    dual_tree_traversal(run_information, fast_sum_interactions, fast_sum_icos_panels);

    if (ID == 0) {
      end = std::chrono::steady_clock::now();
      std::cout << "dual tree traversal time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()
                << " microseconds" << std::endl;
    }
  }



  // compute SAL

  double lon_t, lat_t;

  if (ID == 0) {
    begin = std::chrono::steady_clock::now();
  }

  if (run_information.use_fast) {
    // fast summation
    fast_sum_compute(run_information, fast_sum_interactions, fast_sum_icos_panels, lats, lons, xval, yval, zval, sshs, area, sals, llns);

  } else {
    // direct summation
    for (int i = 0; i < run_information.point_count; i++) {
      if (i % P == ID) {
        lat_t = lats[i];
        lon_t = lons[i];
        for (int j = 0; j < run_information.point_count; j++) {
          sals[i] += sshs[j] * area[j] * sal_ces_gfunc(lat_t, lon_t, lats[j], lons[j], run_information.sph_harm_comps, llns);
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  sync_updates<double>(sals, P, ID, &win_sal, MPI_DOUBLE);

  MPI_Barrier(MPI_COMM_WORLD);

  if (ID == 0) {
    end = std::chrono::steady_clock::now();
    std::cout << "convolution time: " << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()
              << " microseconds" << std::endl;
  }

  std::cout << fast_sum_icos_panels.size() << std::endl;
  std::cout << fast_sum_interactions.size() << std::endl;

  // std::cout << sal_gfunc(0, 0, M_PI / 2.0, 0, 5, llns) << std::endl;

  // write output
  std::string outpath = run_information.out_path + "/" + output_folder + "/output_sal.csv";
  write_state(run_information, sals, outpath);

  MPI_Finalize();
  return 0;
}
