#include <vector>
#include <cmath>

double sbb_coeff(const int deg, const int i, const int j) {
  std::vector<double> log_vals(deg + 1, 0);
  double accum = 0;
  for (int k = 1; k < deg + 1; k++) {
    log_vals[k] = log(k);
    accum += log_vals[k];
  }
  for (int k = 1; k < i + 1; k++) {
    accum -= log_vals[k];
  }
  for (int k = 1; k < j + 1; k++) {
    accum -= log_vals[k];
  }
  for (int k = 1; k < (deg-i-j) + 1; k++) {
    accum -= log_vals[k];
  }
  return exp(accum);
}

void fekete_init(std::vector<std::vector<double>> &points, const int degree) {
  double delta_x = 1.0 / degree;
  int index;
  double a, b, c, part;
  for (int i = 0; i < degree + 1; i++) {
    a = 1 - i * delta_x;
    a = 0.5 * (1 + sin(M_PI / 2 * (2 * a - 1)));
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      c = j * delta_x;
      b = i*delta_x - b;
      b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
      c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
      part = a + b + c;
      points[index][0] = a / part;
      points[index][1] = b / part;
      points[index][2] = c / part;
    }
  }
}

void interp_mat_init_sbb(
    std::vector<double> &mat, const std::vector<std::vector<double>> &points,
    const int degree,
    const int point_count) { // sets up matrix to interpolate with fekete points
  // uses spherical bezier bernstein polynomials
  // for example, for deg 2, evaluates s^2, t^2, u^2, st, su, tu at interpolation points
  int index = 0, place;
  double s, t, u, si, tj, uk, comp, val, spart = 1, tpart=1;
  for (int k = 0; k < point_count; k++) {
    s = points[k][0];
    t = points[k][1];
    u = points[k][2];
    index = 0;
    spart = 1;
    tpart = 1;
    for (int i = 0; i < degree + 1; i++) {
      for (int j = 0; j < degree+1-i; j++) {
        val = spart * tpart;
        if (degree - i -j != 0) {
          val *= pow(u, degree-i-j);
        }
        place = point_count * index + k;
        mat[place] = val;
        index++;
        tpart *= t;
      }
      spart *= s;
    }
  }
}

std::vector<double> interp_vals_sbb(const double s, const double t, const double u, const int degree) {
  // returns vector of SBB basis values of s, t, u
  int count = (degree + 1) * (degree + 2) / 2;
  std::vector<double> out_vals (count, 0);
  double val, factor, spart = 1;
  factor = t / u;
  int index = 0;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    val = spart * pow(u, degree - i);
    for (int j = 0; j < degree + 1 - i; j++) {
      out_vals[index] = val;
      index += 1;
      val *= factor;
    }
    spart *= s;
  }
  return out_vals;
}

double interp_eval_sbb(
    const std::vector<double> &alphas, const double s, const double t, const double u,
    const int degree) { // evaluate SBB interpolation polynomial with coefficients
                        // alpha and barycentric point (s, t, u)
  double accum = 0;
  int index = 0;
  double val, factor, spart = 1;
  factor = t / u;
  for (int i = 0; i < degree + 1; i++) { // degree of s
    val = spart * pow(u, degree - i);
    for (int j = 0; j < degree + 1 - i; j++) {
      accum += val * alphas[index];
      index += 1;
      val *= factor;
    }
    spart *= s;
  }
  return accum;
}
