//
// Created by Riccardo Ceccaroni on 31/01/25.
//

#include <iostream>
#include <vector>
#include <string>
#include "edlib.h"

#include "Correction.h"


void print_job(const std::vector<std::pair<std::string, std::string>>& job) {
    for (const auto& entry : job) {
        std::cout << "Query sequence: " << entry.first << std::endl;
        std::cout << "Target sequence: " << entry.second << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}


std::string run_msa(const std::vector<std::string>& sequences) {
   return sequences[0];
}

std::string correction(const std::vector<std::pair<std::string, std::string>>& job) {
    std::string consensus_sequence;

    for (const auto& entry : job) {
        const std::string& query_sequence = entry.first;
        const std::string& target_sequence = entry.second;

        EdlibAlignResult result = edlibAlign(query_sequence.c_str(), query_sequence.size(),
                                             target_sequence.c_str(), target_sequence.size(),
                                             edlibDefaultAlignConfig());

        if (result.status == EDLIB_STATUS_OK) {
            std::string query_aligned, target_aligned;
            int query_idx = 0, target_idx = 0;

            for (int i = 0; i < result.alignmentLength; ++i) {
                if (result.alignment[i] == EDLIB_EDOP_MATCH || result.alignment[i] == EDLIB_EDOP_MISMATCH) {
                    query_aligned += query_sequence[query_idx++];
                    target_aligned += target_sequence[target_idx++];
                } else if (result.alignment[i] == EDLIB_EDOP_INSERT) {
                    query_aligned += query_sequence[query_idx++];
                    target_aligned += '-';
                } else if (result.alignment[i] == EDLIB_EDOP_DELETE) {
                    query_aligned += '-';
                    target_aligned += target_sequence[target_idx++];
                }
            }

            int t_start = 0;
            int len_window = 10;  // Example window size
            std::vector<std::pair<int, int>> windows_borders;
            int w_id = 0;

            for (size_t i = 0; i < target_aligned.length(); ++i) {
                if (target_aligned[i] != '-') {
                    t_start++;
                    if (t_start % len_window == 0) {
                        w_id = t_start / len_window;
                        windows_borders.emplace_back(i, w_id);
                    }
                }
            }

            w_id++;
            windows_borders.emplace_back(query_aligned.length(), w_id);

            // Perform MSA using ABPOA
            std::vector<std::string> sequences_to_align;
            for (const auto& w : windows_borders) {
                sequences_to_align.push_back(query_sequence.substr(w.first, len_window));
            }
            consensus_sequence += run_msa(sequences_to_align);
        }

        edlibFreeAlignResult(result);
    }

    return consensus_sequence;
}

