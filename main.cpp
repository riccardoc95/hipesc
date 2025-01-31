#include "Jobs.h"
#include <iostream>


#include <vector>
#include <string>

void print_job(const std::vector<std::pair<std::string, std::string>>& job) {
    for (const auto& entry : job) {
        std::cout << "Query sequence: " << entry.first << std::endl;
        std::cout << "Target sequence: " << entry.second << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}


int main() {
    Jobs jobs("data");
    auto job_list = jobs.get_job_list();

    for (size_t i = 0; i < job_list.size(); ++i) {
        std::cout << "Processing job " << i + 1 << " of " << job_list.size() << std::endl;
        auto job = jobs.get(job_list[i]);
        print_job(job);
    }


    return 0;
}