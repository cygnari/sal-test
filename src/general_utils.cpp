#include <vector>
#include <string>
#include <cmath>
#include <tuple>

extern "C" { // lapack
extern int dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
extern int dgelsd_(int *, int *, int *, double *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int *);
extern int dgelss_(int *, int *, int *, double *, int *, double *, int *, double *, double *, int *, double *, int *, int *);
extern int ilaenv_(int *, char *, char *, int *, int *, int *, int *);
}

int linear_solve(const std::vector<double> &a_matrix, std::vector<double> &b_vec, int size, int nrhs, const int solver) {
  // solve Ax=b via dgesv or dgelsd or other options
  std::vector<double> interp_mat_copy = a_matrix; // copy matrix, some methods destroy the matrix
  int info;
  if (solver == 1) { // dgesv
    std::vector<int> ipiv (size, 0);
    dgesv_(&size, &nrhs, &*interp_mat_copy.begin(), &size, &*ipiv.begin(), &*b_vec.begin(), &size, &info);
  } else if (solver == 2) { // dgelss
    std::vector<double> sing_vals(size, 0);
    double rcond = -1;
    int lwork = std::max(3 * size + std::max(2*size, nrhs), 1);
    lwork *= 2;
    std::vector<double> work(lwork, 0);
    int rank;
    dgelss_(&size, &size, &nrhs, &*interp_mat_copy.begin(), &size, &*b_vec.begin(),
        &size, &*sing_vals.begin(), &rcond, &rank, &*work.begin(), &lwork, &info);
  } else if (solver == 3) { // dgelsd
    std::vector<double> sing_vals(size, 0);
    int rank;
    double rcond = -1;
    std::string name = "DGELSD";
    std::string empty = " ";
    int nine = 9;
    int zero = 0;
    int smallsize = ilaenv_(&nine, &name[0], &empty[0], &zero, &zero, &zero, &zero);
    int nlvl = std::max(1, int(log2(size/(smallsize+1)))+1);
    int lwork = 12*size+2*size*smallsize+8*size*nlvl+size*nrhs+pow(smallsize+1,2);
    lwork *= 2;
    int liwork = std::max(1, 3*size*nlvl+11*size);
    std::vector<double> work(lwork, 0);
    std::vector<int> iwork(liwork, 0);
    dgelsd_(&size, &size, &nrhs, &*interp_mat_copy.begin(), &size, &*b_vec.begin(), &size, &*sing_vals.begin(), &rcond, &rank, &*work.begin(), &lwork, &*iwork.begin(), &info);
  }
  return info;
}

std::tuple<double, double, double> latlon_to_xyz(const double lat, const double lon, const double radius) {
  // converts latitude, longitude to x y z coords
  std::tuple<double, double, double> xyz;
  double x, y, z, dist;
  x = radius * cos(lat) * cos(lon);
  y = radius * cos(lat) * sin(lon);
  z = radius * sin(lat);
  dist = sqrt(x * x + y * y + z * z);
  xyz = std::make_tuple(x / dist, y / dist, z / dist);
  return xyz;
}

std::vector<double> project_to_sphere(double x, double y, double z, const double radius) {
  // projects (x, y, z) to sphere of radius
  double point_dist = sqrt(x * x + y * y + z * z);
  x /= point_dist;
  y /= point_dist;
  z /= point_dist;
  return std::vector<double> {x, y, z};
}
