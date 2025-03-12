//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#include "correction.h"
#include "consent.h"


std::string correction(std::string& query, std::vector<Overlap>& targets, const std::string& method) {
    if (method == "default") {
        return consent_correction(query, targets);
    }else {
        return "";
    }
}
