#include <vector>
#include "structs.hpp"
#include "green_func.hpp"
#include "general_utils.hpp"
#include "interp_utils.hpp"

void pp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns) {
  // compute particle particle interaction
  int target_i, source_j, source_count;
  double ssh, lon_x, lon_y, lat_x, lat_y, gfval;

  source_count = icos_panels[index_source].point_count;

  std::vector<int> sources_j (source_count, 0);
  std::vector<double> lons_y (source_count, 0), lats_y (source_count, 0), area_y (source_count, 0), sshs_y (source_count, 0);

  for (int j = 0; j < source_count; j++) {
    source_j = icos_panels[index_source].points_inside[j];
    sources_j[j] = source_j;
    lons_y[j] = lons[source_j];
    lats_y[j] = lats[source_j];
    area_y[j] = area[source_j];
    sshs_y[j] = sshs[source_j];
  }

  for (int i = 0; i < icos_panels[index_target].point_count; i++) {
    target_i = icos_panels[index_target].points_inside[i];
    lon_x = lons[target_i];
    lat_x = lats[target_i];
    for (int j = 0; j < source_count; j++) {
      gfval = sal_gfunc(lat_x, lon_x, lats_y[j], lons_y[j], run_information.sph_harm_comps, llns);
      sals[target_i] += gfval * area_y[j] * sshs_y[j];
    }
  }
}

void pc_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns) {
  // computes particle cluster interaction
  int iv1s, iv2s, iv3s, point_index, offset, info, dim = run_information.interp_point_count, count_target = icos_panels[index_target].point_count, count_source = icos_panels[index_source].point_count;
  std::vector<double> v1s, v2s, v3s, target_particle, bary_cord, source_particle, interp_matrix(dim * dim, 0), proxy_weights(dim, 0),
                      basis_vals, func_vals(dim * count_target, 0), tlons (count_target, 0), tlats (count_target, 0);
  double ssh, us, vs, ws, tlat, tlon, cx, cy, cz, clat, clon, sx, sy, sz, alpha, scalar;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  fekete_init(interp_points, run_information.interp_degree);
  v1s = icos_panels[index_source].vertex_1;
  v2s = icos_panels[index_source].vertex_2;
  v3s = icos_panels[index_source].vertex_3;

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    tlons[i] = lons[point_index];
    tlats[i] = lats[point_index];
  }

  for (int i = 0; i < count_source; i++) { // compute proxy weights
    point_index = icos_panels[index_source].points_inside[i];
    sx = xs[point_index];
    sy = ys[point_index];
    sz = zs[point_index];
    ssh = sshs[point_index];
    bary_cord = barycoords(v1s, v2s, v3s, sx, sy, sz);
    basis_vals = interp_vals_sbb(bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
    for (int j = 0; j < dim; j++) {
      proxy_weights[j] += basis_vals[j] * ssh * area[point_index];
    }
  }

  for (int i = 0; i < dim; i++) { // set up interpolation matrix
    us = interp_points[i][0];
    vs = interp_points[i][1];
    ws = 1.0 - us - vs;
    cx = us * v1s[0] + vs * v2s[0] + ws * v3s[0];
    cy = us * v1s[1] + vs * v2s[1] + ws * v3s[1];
    cz = us * v1s[2] + vs * v2s[2] + ws * v3s[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1s, v2s, v3s, cx, cy, cz);
    interp_points[i]=bary_cord;
    auto [clat, clon] = xyz_to_latlon(cx, cy, cz);
    for (int j = 0; j < count_target; j++) {
      tlon = tlons[j];
      tlat = tlats[j];
      offset = j * dim + i;
      func_vals[offset] = sal_gfunc(tlat, tlon, clat, clon, run_information.sph_harm_comps, llns);
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  info = linear_solve(interp_matrix, func_vals, dim, count_target, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in pc vel computation");
  }

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    for (int j = 0; j < dim; j++) {
      alpha = func_vals[dim*i+j];
      sals[point_index] += alpha * proxy_weights[j];
    }
  }
}

void cp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns) {
  // computes cluster particle interaction
  int iv1, iv2, iv3, point_index, dim = run_information.interp_point_count, count_target = icos_panels[index_target].point_count, count_source = icos_panels[index_source].point_count;
  std::vector<double> v1, v2, v3, bary_cord, interptargets(dim, 0), interp_matrix(dim * dim, 0), lons_y (count_source, 0), lats_y (count_source, 0), area_y (count_source, 0), sshs_y (count_source, 0);
  double u, v, w, vor, cx, cy, cz, slon, slat, clon, clat, tx, ty, tz, scalar;
  std::vector<std::vector<double>> interp_points(dim, std::vector<double>(3, 0));

  for (int j = 0; j < count_source; j++) {
    point_index = icos_panels[index_source].points_inside[j];
    lons_y[j] = lons[point_index];
    lats_y[j] = lats[point_index];
    area_y[j] = area[point_index];
    sshs_y[j] = sshs[point_index];
  }

  fekete_init(interp_points, run_information.interp_degree);

  v1 = icos_panels[index_target].vertex_1;
  v2 = icos_panels[index_target].vertex_2;
  v3 = icos_panels[index_target].vertex_3;
  for (int i = 0; i < dim; i++) {
    u = interp_points[i][0];
    v = interp_points[i][1];
    w = 1.0 - u - v;
    cx = u * v1[0] + v * v2[0] + w * v3[0];
    cy = u * v1[1] + v * v2[1] + w * v3[1];
    cz = u * v1[2] + v * v2[2] + w * v3[2];
    scalar = run_information.radius / sqrt(cx * cx + cy * cy + cz * cz);
    cx *= scalar;
    cy *= scalar;
    cz *= scalar;
    bary_cord = barycoords(v1, v2, v3, cx, cy, cz);
    interp_points[i] = bary_cord;
    auto [clat, clon] = xyz_to_latlon(cx, cy, cz);
    for (int j = 0; j < count_source; j++) {
      point_index = icos_panels[index_source].points_inside[j];
      slon = lons_y[j];
      slat = lats_y[j];
      interptargets[i] += sal_gfunc(clat, clon, slat, slon, run_information.sph_harm_comps, llns) * sshs_y[j] * area_y[j];
    }
  }

  interp_mat_init_sbb(interp_matrix, interp_points, run_information.interp_degree, dim);

  int info = linear_solve(interp_matrix, interptargets, dim, 1, 3);
  if (info > 0) {
    throw std::runtime_error("Error with linear solve in cp vel computation");
  }

  for (int i = 0; i < count_target; i++) {
    point_index = icos_panels[index_target].points_inside[i];
    tx = xs[point_index];
    ty = ys[point_index];
    tz = zs[point_index];
    bary_cord = barycoords(v1, v2, v3, tx, ty, tz);
    sals[point_index] += interp_eval_sbb(interptargets, bary_cord[0], bary_cord[1], bary_cord[2], run_information.interp_degree);
  }
}

void cc_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns) {
  // computes cluster cluster interaction
}

void fast_sum_compute(const RunConfig& run_information, const std::vector<InteractPair> interactions, const std::vector<IcosPanel>& icos_panels,
                      const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                      const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns) {
  // computes all the interactions
  for (int i = 0; i < interactions.size(); i++) {
    if (interactions[i].interact_type == 0) {
      // PP
      pp_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, xs, ys, zs, sshs, area, sals, llns);
    } else if (interactions[i].interact_type == 1) {
      // PC
      pc_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, xs, ys, zs, sshs, area, sals, llns);
    } else if (interactions[i].interact_type == 2) {
      // CP
      cp_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, xs, ys, zs, sshs, area, sals, llns);
    } else {
      // CC
      cc_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, xs, ys, zs, sshs, area, sals, llns);
    }
  }
}
