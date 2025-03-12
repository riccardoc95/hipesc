//
// Created by Riccardo Ceccaroni on 12/03/25.
//

#ifndef UTILS_H
#define UTILS_H

struct Overlap {
    std::string target;
    size_t query_start;
    size_t query_end;
    size_t target_start;
    size_t target_end;
    std::string strand;
};

#endif //UTILS_H
