//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <algorithm>
#include <sys/resource.h>

struct Target {
    std::string target_name;
    size_t query_start;
    size_t query_end;
    size_t target_start;
    size_t target_end;
    std::string strand;
};

struct Job {
    std::string query_name;
    std::vector<Target> jobs;
};

struct Alignment {
    int query_start;
    int query_end;
    std::string target_name;
    int target_start;
    int target_end;
};

void print_memory_usage();
std::string datetime();
void trim_newline(std::string &str);
void print_time_message(const std::string& message, int rank=0);

#endif //UTILS_H
