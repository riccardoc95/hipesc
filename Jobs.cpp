//
// Created by Riccardo Ceccaroni on 30/01/25.
//
#include <iostream>
#include "Jobs.h"

// Jobs class (new part to integrate FastQ and Paf)

Jobs::Jobs(const std::string& folder) : fastq(folder), paf(folder) {}

std::vector<std::string> Jobs::get_job_list() {
    std::vector<std::string> job_list;
    for (const auto& pair : paf.index) {
        job_list.push_back(pair.first);
    }
    return job_list;
}

std::vector<std::pair<std::string, std::string>> Jobs::get(const std::string& query_name) {
    std::vector<std::pair<std::string, std::string>> job;
    auto paf_results = paf.get(query_name);
    for (const auto& entry : paf_results) {
        const auto& query_name = entry.first;
        auto query = fastq.get(query_name);
        for (const auto& tuple : entry.second) {
            const auto& target_name = std::get<0>(tuple);
            const auto& query_start = std::get<1>(tuple);
            const auto& query_end = std::get<2>(tuple);
            const auto& target_start = std::get<3>(tuple);
            const auto& target_end = std::get<4>(tuple);
            auto target = fastq.get(target_name);

            // Corrected substring extraction: calculate length
            int q_start = std::stoi(query_start);
            int q_end = std::stoi(query_end);
            int t_start = std::stoi(target_start);
            int t_end = std::stoi(target_end);

            job.push_back({query["sequence"].substr(q_start, q_end - q_start),
                           target["sequence"].substr(t_start, t_end - t_start)});
        }
    }
    return job;
}



