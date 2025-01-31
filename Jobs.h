//
// Created by Riccardo Ceccaroni on 30/01/25.
//

#ifndef JOBS_H
#define JOBS_H

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include "FastQ.h"
#include "Paf.h"

class Jobs {
public:
    explicit Jobs(const std::string& folder);
    std::vector<std::string> get_job_list();
    std::vector<std::pair<std::string, std::string>> get(const std::string& query_name);

//private:
    FastQ fastq;
    Paf paf;
};


#endif //JOBS_H
