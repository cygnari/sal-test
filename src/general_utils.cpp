#include <vector>
#include <string>
#include <cmath>
#include <tuple>
#include <cassert>

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

std::tuple<double, double> xyz_to_latlon(const double x, const double y, const double z) {
  // turns cartesian coordinates to spherical coordinates
  double colat, lon;
  colat = atan2(sqrt(x * x + y * y), z); // colatitude
  lon = atan2(y, x);                   // longitude
  return std::make_tuple(M_PI / 2.0 - colat, lon);
}

std::tuple<double, double> xyz_to_latlon(const std::vector<double>& point) {
  // turns cartesian coordinates to spherical coordinates
  double colat, lon;
  colat = atan2(sqrt(point[0] * point[0] + point[1] * point[1]), point[2]); // colatitude
  lon = atan2(point[1], point[2]);                                          // longitude
  return std::make_tuple(M_PI / 2.0 - colat, lon);
}

std::vector<double> project_to_sphere(double x, double y, double z, const double radius) {
  // projects (x, y, z) to sphere of radius
  double point_dist = sqrt(x * x + y * y + z * z);
  x /= point_dist;
  y /= point_dist;
  z /= point_dist;
  return std::vector<double> {x, y, z};
}

double gcdist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double radius) {
  double s = (x1 * x2 + y1 * y2 + z1 * z2) / (sqrt(x1 * x1 + y1 * y1 + z1 * z1) * sqrt(x2 * x2 + y2 * y2 + z2 * z2));
  double theta = acos(std::min(std::max(s, -1.0), 1.0));
  return theta * radius;
}

double gcdist(const double lat1, const double lon1, const double lat2, const double lon2, const double radius) {
  return radius * acos(std::min(1.0, std::max(-1.0, sin(lat1) * sin(lat2) +
                                cos(lat1) * cos(lat2) * cos(lon2 - lon1))));
}

std::vector<double> barycoords(const std::vector<double> &p1,
                               const std::vector<double> &p2,
                               const std::vector<double> &p3,
                               const double x, const double y, const double z) {
  // finds triangle barycentric coordinates of point p
  assert(p1.size() == 3);
  assert(p2.size() == 3);
  assert(p3.size() == 3);
  std::vector<double> coords {x, y, z};
  std::vector<double> mat{p1[0], p1[1], p1[2], p2[0], p2[1],
                          p2[2], p3[0], p3[1], p3[2]};
  int info = linear_solve(mat, coords, 3, 1, 1);
  if (info != 0) {
    throw std::runtime_error("Error in barycentric coordinate computation, line 94");
  }
  return coords;
}
