#ifndef H_GREEN_FUNC_H
#define H_GREEN_FUNC_H

double sal_gfunc(const double lat_t, const double lon_t, const double lat_s, const double lon_s,
                 const int sph_comps, const std::vector<double>& llns);

double sal_ces_gfunc(const double lat_t, const double lon_t, const double lat_s, const double lon_s,
                 const int sph_comps, const std::vector<double>& llns);

#endif
