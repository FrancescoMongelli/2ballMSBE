#ifndef LSH_H
#define LSH_H

#include <vector>
#include <unordered_set>
#include "../Hash.cpp"

struct Biclique {
    int id;                 // global id
    bool type;              // 0 = 1 hop bicliques, 1 = 2 hop bicliques
    vector<uint32_t> nodes; // biclique nodes
};

struct pair_hash {
    size_t operator()(const pair<int,int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

uint32_t* computeMinHashSignature(
    const Biclique& b,
    std::vector<Hash<uint32_t>*>& hashes
);

std::unordered_set<std::pair<int,int>, pair_hash>
runLSH(
    const std::vector<Biclique>& bicliques,
    std::vector<Hash<uint32_t>*>& hashes,
    int r,
    int b
);

double dissimilarityIndex(
    const std::unordered_set<std::pair<int,int>, pair_hash>& pairs,
    int nA,
    int nB
);

void avgDissimilarity(
    const string &dispath,
    ofstream &out
);

#endif
