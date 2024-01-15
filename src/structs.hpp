#ifndef H_STRUCTS_H
#define H_STRUCTS_H

#include <string>
#include <vector>

struct RunConfig {
  bool use_fast = false;
  bool use_cesaro = false; // use cesaro summed kernel or not
  std::string out_path;    // ../run-output/ locally, on Derecho, /glade/derecho/scratch/achen/bve/
  int write_precision = 6; // number of decimal places, 6 for data visualization, 16 for error testing
  bool write_output = false;
  double radius = 1.0;
  int point_count;    // number of dynamics points, full point count is 2313486
  int interp_degree;      // interpolation degree
  int interp_point_count; // number of interpolation points
  int info_per_point; // how many doubles each point is, for example, storing x y z vor tracer = 5
  int time; // time 0, 1, 2, or 3
  int sph_harm_comps; // number of spherical harmonic components to use

  // fast sum info
  int fast_sum_cluster_thresh; // threshold for a triangle being a cluster
  double fast_sum_theta;       // well separated threshold
  bool fast_sum_rotate = true; // whether or not to rotate the fast sum grid
  double fast_sum_rotate_alph = 0.01; // 3 rotation coefficients
  double fast_sum_rotate_beta = 0.01;
  double fast_sum_rotate_gamm = 0.01;

  // mpi info
  int mpi_P;       // total MPI ranks
  int mpi_ID;      // own MPI rank
  int particle_lb; // range of assigned particles
  int particle_ub;
  int particle_own;
  int target_lb;
  int target_ub;
  int target_own;
};

struct IcosPanel {
  int level = 0; // refinement level, level 0 = base icosahedron
  bool is_leaf = true; // whether or not this is a leaf panel
  int id; // location in vector
  IcosPanel* parent_panel;
  IcosPanel* child_panel_1;
  IcosPanel* child_panel_2;
  IcosPanel* child_panel_3;
  IcosPanel* child_panel_4;
  std::vector<double> vertex_1 {0, 0, 0};
  std::vector<double> vertex_2 {0, 0, 0};
  std::vector<double> vertex_3 {0, 0, 0};
  std::vector<double> center_p {0, 0, 0};
  double radius;
  int point_count = 0; // points inside the panel
  std::vector<int> points_inside;
};

struct InteractPair {
  int index_target;
  int index_source;
  int interact_type; // 0 for PP, 1 for PC, 2 for CP, 3 for CC
};

#endif
