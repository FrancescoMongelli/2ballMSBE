#ifndef LSH_CPP
#define LSH_CPP

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "./include/LSH.h"

using namespace std;

uint32_t* computeMinHashSignature(
    const Biclique& b,
    vector<Hash<uint32_t>*>& hashes
) {
    int K = hashes.size();
    uint32_t* sig = new uint32_t[K];

    for (int i = 0; i < K; i++) {
        uint32_t minv = UINT32_MAX;
        for (uint32_t x : b.nodes)
            minv = min(minv, (*hashes[i])(x));
        sig[i] = minv;
    }
    return sig;
}

string bandKey(uint32_t* sig, int start, int r) {
    ostringstream oss;
    for (int i = 0; i < r; i++) {
        oss << sig[start + i];
        if (i + 1 < r) oss << ",";
    }
    return oss.str();
}

vector<unordered_multimap<string,int>>
buildLSHBuckets(uint32_t** signatures, int n, int r, int b) {

    vector<unordered_multimap<string,int>> buckets(b);

    for (int i = 0; i < n; i++) {
        for (int band = 0; band < b; band++) {
            string key = bandKey(signatures[i], band*r, r);
            buckets[band].insert({key, i});
        }
    }
    return buckets;
}

unordered_set<pair<int,int>, pair_hash>
extractCandidatePairs(
    const vector<unordered_multimap<string,int>>& buckets,
    const vector<Biclique>& bicliques
) {
    unordered_set<pair<int,int>, pair_hash> candidates;

    for (const auto& band : buckets) {

        unordered_map<string, vector<int>> localBuckets;

        for (auto& it : band)
            localBuckets[it.first].push_back(it.second);

        for (auto& [key, ids] : localBuckets) {
            if (ids.size() < 2) continue;

            for (size_t i = 0; i < ids.size(); i++) {
                for (size_t j = i + 1; j < ids.size(); j++) {

                    int a = ids[i];
                    int b = ids[j];

                    // only cross-file
                    if (bicliques[a].type != bicliques[b].type) {
                        if (a > b) swap(a, b);
                        candidates.insert({a, b});
                    }
                }
            }
        }
    }
    return candidates;
}

unordered_set<pair<int,int>, pair_hash>
runLSH(
    const vector<Biclique>& bicliques,
    vector<Hash<uint32_t>*>& hashes,
    int r,
    int b
) {
    int n = bicliques.size();

    uint32_t** signatures = new uint32_t*[n];
    for (int i = 0; i < n; i++)
        signatures[i] = computeMinHashSignature(bicliques[i], hashes);

    auto buckets = buildLSHBuckets(signatures, n, r, b);
    auto candidates = extractCandidatePairs(buckets, bicliques);

    for (int i = 0; i < n; i++)
        delete[] signatures[i];
    delete[] signatures;

    return candidates;
}

double dissimilarityIndex(
    const unordered_set<pair<int,int>, pair_hash>& candidates,
    int nA,
    int nB
) {
    if (nA == 0 || nB == 0) return 1.0;

    vector<bool> coveredA(nA, false);
    vector<bool> coveredB(nB, false);

    for (const auto& [i, j] : candidates) {
        // i ∈ A, j ∈ B
        coveredA[i] = true;
        coveredB[j - nA] = true;
    }

    int covA = 0, covB = 0;
    for (bool x : coveredA) if (x) covA++;
    for (bool x : coveredB) if (x) covB++;

    double similarity =
        double(covA + covB) / double(nA + nB);

    return 1.0 - similarity;
}


void avgDissimilarity(
    const string &dispath,
    ofstream &out
) {
    ifstream in(dispath);

    if (!in.is_open()) {
        cerr << "Error while opening file: " << dispath << endl << flush;
        return;
    }

    string header;
    getline(in, header);

    int run = 0;
    std::vector<double> dissimilarities;
    double dissim;
    double sum = 0.0;

    while (in >> run >> dissim) {
        dissimilarities.push_back(dissim);
    }

    in.close();

    if (run == 0) {
        cerr << "Nessun dato valido in " << dispath << endl;
        return;
    }

    for (double dis: dissimilarities) {
        sum += dis;
    }

    double avg = sum / run;

    double std_dev = 0.0;

    for (double dis: dissimilarities) {
        std_dev += pow(avg - dis, 2);
    }

    std_dev = sqrt(std_dev/run); 

    out << run << " " << avg << " " << std_dev << endl;
}

#endif