#ifndef H_IO_UTILS_H
#define H_IO_UTILS_H

#include "structs.hpp"
#include <string>

void read_run_config(const std::string file_name, RunConfig &run_information);

std::string create_config(const RunConfig &run_information);

void read_data(const RunConfig &run_information, std::vector<double>& latvals, std::vector<double>& lonvals, std::vector<double>& sshs,
               std::vector<double>& areas);

void load_llns(const RunConfig &run_information, std::vector<double>& llns);

void write_state(const RunConfig &run_information, const std::vector<double> &sal, const std::string path);

#endif
