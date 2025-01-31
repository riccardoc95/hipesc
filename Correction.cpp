//
// Created by Riccardo Ceccaroni on 31/01/25.
//

#include "Correction.h"


void print_job(const std::vector<std::pair<std::string, std::string>>& job) {
    for (const auto& entry : job) {
        std::cout << "Query sequence: " << entry.first << std::endl;
        std::cout << "Target sequence: " << entry.second << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}


void align_and_process_windows(const std::vector<std::pair<std::string, std::string>>& job) {
    for (const auto& entry : job) {
        std::cout << "Query sequence: " << entry.first << std::endl;
        std::cout << "Target sequence: " << entry.second << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}