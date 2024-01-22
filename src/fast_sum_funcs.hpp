#ifndef H_FAST_SUM_FUNCS_H
#define H_FAST_SUM_FUNCS_H

#include <vector>
#include "structs.hpp"

void pp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns);

void pc_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns);

void cp_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns);

void cc_interaction(const RunConfig& run_information, const int index_target, const int index_source, const std::vector<IcosPanel>& icos_panels,
                    const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                    const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns);

void fast_sum_compute(const RunConfig& run_information, const std::vector<InteractPair> interactions, const std::vector<IcosPanel>& icos_panels,
                      const std::vector<double>& lats, const std::vector<double>& lons, const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs,
                      const std::vector<double>& sshs, const std::vector<double>& area, std::vector<double>& sals, const std::vector<double>& llns);

#endif
