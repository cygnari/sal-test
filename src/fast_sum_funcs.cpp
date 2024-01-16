#include <vector>
#include "structs.hpp"
#include "green_func.hpp"

void pp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& sshs, const std::vector<double>& area,
                    std::vector<double>& sals, const std::vector<double>& llns) {
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
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& sshs, const std::vector<double>& area,
                    std::vector<double>& sals, const std::vector<double>& llns) {
  // computes particle cluster interaction
}

void cp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& sshs, const std::vector<double>& area,
                    std::vector<double>& sals, const std::vector<double>& llns) {
  // computes cluster particle interaction
}

void cc_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& sshs, const std::vector<double>& area,
                    std::vector<double>& sals, const std::vector<double>& llns) {
  // computes cluster cluster interaction
}

void fast_sum_compute(const RunConfig& run_information, const std::vector<InteractPair> interactions, const std::vector<IcosPanel>& icos_panels,
                      const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& sshs, const std::vector<double>& area,
                      std::vector<double>& sals, const std::vector<double>& llns) {
  // computes all the interactions
  for (int i = 0; i < interactions.size(); i++) {
    if (interactions[i].interact_type == 0) {
      // PP
      pp_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, sshs, area, sals, llns);
    } else if (interactions[i].interact_type == 1) {
      // PC
      pc_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, sshs, area, sals, llns);
    } else if (interactions[i].interact_type == 2) {
      // CP
      cp_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, sshs, area, sals, llns);
    } else {
      // CC
      cc_interaction(run_information, interactions[i].index_target, interactions[i].index_source, icos_panels, lats, lons, sshs, area, sals, llns);
    }
  }
}
