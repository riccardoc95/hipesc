//
// Created by Riccardo Ceccaroni on 31/01/25.
//

#ifndef CORRECTION_H
#define CORRECTION_H

#include <vector>
#include <string>
#include <iostream>

// Function declaration
void print_job(const std::vector<std::pair<std::string, std::string>>& job);
std::string run_msa(const std::vector<std::string>& sequences);
std::string correction(const std::vector<std::pair<std::string, std::string>>& job);




#endif //CORRECTION_H
