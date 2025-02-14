//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <string>

std::string compress_zstd(const std::string& str);
std::string decompress_zstd(const std::string& compressedStr, size_t originalSize);

#endif //COMPRESSION_H
