//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#ifndef PAF_H
#define PAF_H

#include <mpi.h>
#include <omp.h>
#include <string>

void parse_paf_line(std::string line_str, std::string& query_name, std::string& target_name, size_t& query_start, size_t& query_end, size_t& target_start, size_t& target_end, std::string& strand);

#endif //PAF_H
