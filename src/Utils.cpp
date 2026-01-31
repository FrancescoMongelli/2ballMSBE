#ifndef UTILS_CPP
#define UTILS_CPP

#include "include/Utils.h"
#include "include/LSH.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>

using namespace std;

uint32_t *read_edges(string fname, uint32_t *n1, uint32_t* n2, uint64_t *m)
{
    ifstream file(fname);

    file >> *n1;
    file >> *n2;
    file >> *m;

    uint32_t *edges = new uint32_t[2 * (*m)];
    uint32_t from, to;
    uint64_t i = 0;

    while (file >> from)
    {
        file >> to;
        edges[i++] = from;
        edges[i++] = to;
    }

    file.clear();
    file.close();

    return edges;
}

vector<uint32_t> parseLine(const string& line)
{
    vector<uint32_t> v;
    stringstream ss(line);
    uint32_t x;
    while (ss >> x)
        v.push_back(x);
    return v;
}

vector<Biclique> parseBicliqueFile(const string& path, bool typeId)
{
    ifstream file(path);
    vector<Biclique> bicliques;

    if (!file.is_open())
        throw runtime_error("Cannot open file: " + path);

    string line;
    int id = 0;

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        if (line.rfind("eps:", 0) == 0 ||
            line.rfind("tau:", 0) == 0 ||
            line.rfind("k", 0) == 0 ||
            line.rfind("phi", 0) == 0)
            continue;

        // left nodes
        vector<uint32_t> left = parseLine(line);

        // right nodes
        string line2;
        getline(file, line2);
        vector<uint32_t> right = parseLine(line2);

        // full biclique
        unordered_set<uint32_t> nodeset;
        for (auto x : left)  nodeset.insert(x);
        for (auto x : right) nodeset.insert(x);

        Biclique b;
        b.id = id++;
        b.type = typeId;
        b.nodes.assign(nodeset.begin(), nodeset.end());

        bicliques.push_back(std::move(b));
    }

    return bicliques;
}


#endif