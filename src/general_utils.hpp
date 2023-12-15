#ifndef H_GENERAL_UTIL_H
#define H_GENERAL_UTIL_H

#include <vector>

int linear_solve(const std::vector<double> &a_matrix, std::vector<double> &b_vec, int size, int nrhs, const int solver);

std::tuple<double, double, double> latlon_to_xyz(const double lat, const double lon, const double radius);

std::vector<double> project_to_sphere(double x, double y, double z, const double radius);

#endif
