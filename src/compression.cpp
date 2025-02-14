//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#include "compression.h"

#include <zstd.h>
#include <vector>
#include <stdexcept>

std::string compress_zstd(const std::string& str) {
    size_t compressedSize = ZSTD_compressBound(str.size());
    std::vector<char> compressedData(compressedSize);

    compressedSize = ZSTD_compress(compressedData.data(), compressedData.size(), str.data(), str.size(), 3);

    if (ZSTD_isError(compressedSize)) {
        throw std::runtime_error("Zstd Compression failed");
    }

    return compressedData.data();
}

std::string decompress_zstd(const std::string& compressedStr, size_t originalSize) {
    std::vector<char> decompressedData(originalSize);

    size_t decompressedSize = ZSTD_decompress(decompressedData.data(), decompressedData.size(), compressedStr.data(), compressedStr.size());

    if (ZSTD_isError(decompressedSize)) {
        throw std::runtime_error("Zstd Decompression failed");
    }

    return decompressedData.data();
}
