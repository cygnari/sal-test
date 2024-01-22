#ifndef H_INTERP_UTILS_H
#define H_INTERP_UTILS_H

#include "structs.hpp"

double sbb_coeff(const int deg, const int i, const int j);

void fekete_init(std::vector<std::vector<double>> &points, const int degree);

void interp_mat_init_sbb(std::vector<double> &mat, const std::vector<std::vector<double>> &points,
   const int degree, const int point_count);

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree);

double interp_eval_sbb(const std::vector<double> &alphas, const double s, const double t, const double u,
   const int degree);

#endif
