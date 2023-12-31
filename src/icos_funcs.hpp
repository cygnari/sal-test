#ifndef H_ICOS_FUNCS_H
#define H_ICOS_FUNCS_H

#include <vector>
#include "structs.hpp"

void initialize_icosahedron(const RunConfig& run_information, std::vector<IcosPanel>& icos_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z);

#endif
