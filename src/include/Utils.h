#ifndef UTILS_H
#define UTILS_H

#include <cstdint>
#include <string>
#include <vector>

#include "LSH.h"

uint32_t *read_edges(std::string fname, uint32_t *n1, uint32_t *n2, uint64_t *m);
std::vector<uint32_t> parseLine(const std::string& line);
std::vector<Biclique> parseBicliqueFile(const std::string& path, bool typeId);

#endif