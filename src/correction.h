//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#ifndef CORRECTION_H
#define CORRECTION_H

#include <string>
#include <vector>

struct Overlap {
    std::string target;
    size_t query_start;
    size_t query_end;
    size_t target_start;
    size_t target_end;
    std::string strand;
};

std::string correction(std::string& query, std::vector<Overlap>& targets, const std::string& method);

#endif //CORRECTION_H
