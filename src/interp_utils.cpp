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
    part = a;
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      c = j * delta_x;
      b = i*delta_x - b;
      b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
      c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
      part += b + c;
      points[index][0] = a / part;
      points[index][1] = b / part;
      points[index][2] = c / part;
    }
  }
}

void interp_mat_init(
    std::vector<double> &mat, const std::vector<std::vector<double>> &points,
    const int degree,
    const int point_count) { // sets up matrix to interpolate with fekete points
  // simple polynomial in s and t, barycentric coordinates
  // for example, for deg 2, evaluates 1, s, t, s^2, st, t^2 at interpolation points
  int index, place;
  double a, b;
  for (int i = 0; i < degree + 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      for (int k = 0; k < point_count; k++) {
        a = points[k][0];
        b = points[k][1];
        place = index * point_count + k;
        mat[place] = pow(a, i - j) * pow(b, j);
      }
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
  double s, t, u, si, tj, uk, comp, val, spart = 1, factor, tpart=1;
  for (int k = 0; k < point_count; k++) {
    s = points[k][0];
    t = points[k][1];
    u = points[k][2];
    index = 0;
    spart = 1;
    // factor = t/u;
    for (int i = 0; i < degree + 1; i++) {
      // val = spart * pow(u, degree-i);
      for (int j = 0; j < degree+1-i; j++) {
        val = spart * tpart * pow(u, degree-i-j);
        place = point_count * index + k;
        mat[place] = val;
        index++;
        val *= factor;
        tpart *= t;
      }
      spart *= s;
    }
  }
}

double interp_eval(
    const std::vector<double> &alphas, const double s, const double t,
    const int degree) { // evaluate interpolation polynomial with coefficients
                        // alpha and barycentric point (s, t)
  double accum = 0;
  int index;
  for (int i = 0; i < degree + 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      index = i * (i + 1) / 2 + j;
      accum += pow(s, i - j) * pow(t, j) * alphas[index];
    }
  }
  return accum;
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
