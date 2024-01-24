#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

double sal_gfunc(const double lat_t, const double lon_t, const double lat_s, const double lon_s,
                 const int sph_comps, const std::vector<double>& llns) {
  double coeff = 3 * 1.035 / (5.517 * sqrt(4 * M_PI));
  // first compute the angle between source and target
  double cos_angle = std::min(1.0, std::max(-1.0, sin(lat_t) * sin(lat_s) +
                                cos(lat_t) * cos(lat_s) * cos(lon_t - lon_s)));

  // next compute all the needed spherical harmonic components
  std::vector<double> legendre_poly_vals (sph_comps, 0);
  if (sph_comps >= 1) {
    legendre_poly_vals[0] = 1;
  }
  if (sph_comps >= 2) {
    legendre_poly_vals[1] = cos_angle;
  }
  if (sph_comps >= 3) {
    for (int i = 2; i < sph_comps; i++) {
      legendre_poly_vals[i] = ((2 * i - 1) * cos_angle * legendre_poly_vals[i-1] - (i - 1) * legendre_poly_vals[i-2]) / (1.0 * i);
    }
  }

  // for (int i = 0; i < sph_comps; i++) {
  //   std::cout << legendre_poly_vals[i] * llns[i] * coeff / sqrt(2 * i + 1) << std::endl;
  // }


  double output = 0;
  for (int i = 0; i < sph_comps; i++) {
    output += legendre_poly_vals[i] * llns[i] / sqrt(2 * i + 1);
  }
  return output * coeff;
}

double sal_ces_gfunc(const double lat_t, const double lon_t, const double lat_s, const double lon_s,
                 const int sph_comps, const std::vector<double>& llns) {
  // first compute the angle between source and target
  double cos_angle = std::min(1.0, std::max(-1.0, sin(lat_t) * sin(lat_s) +
                                cos(lat_t) * cos(lat_s) * cos(lon_t - lon_s)));

  // next compute all the needed spherical harmonic components
  std::vector<double> legendre_poly_vals (sph_comps, 0);
  if (sph_comps >= 1) {
    legendre_poly_vals[0] = 1;
  }
  if (sph_comps >= 2) {
    legendre_poly_vals[1] = cos_angle;
  }
  if (sph_comps >= 3) {
    for (int i = 2; i < sph_comps; i++) {
      legendre_poly_vals[i] = ((2 * i - 1) * cos_angle * legendre_poly_vals[i-1] - (i - 1) * legendre_poly_vals[i-2]) / (1.0 * i);
    }
  }

  double coeff = 3 * 1.035 / (5.517 * sqrt(4 * M_PI));
  double output = 0;
  for (int i = 0; i < sph_comps; i++) {
    output += legendre_poly_vals[i] * llns[i] / sqrt(2 * i + 1) * (1 - i / sph_comps);
  }
  return output * coeff;
}
