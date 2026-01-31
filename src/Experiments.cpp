#ifndef EXPERIMENTS_CPP
#define EXPERIMENTS_CPP

#include "include/Experiments.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int numberOfBicliques(string fpath)
{
    ifstream file(fpath);

    string line;
    int validLines = 0;

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        if (line.rfind("eps:", 0) == 0 || line.rfind("tau:", 0) == 0 || line.rfind("k", 0) == 0 || line.rfind("phi", 0) == 0)
            continue;

        validLines++;
    }

    int bicliqueCount = validLines / 2;

    return bicliqueCount;
}

int countSize(string line)
{
    stringstream ss(line);
    int x, count = 0;
    while (ss >> x)
        count++;
    return count;
}

int maxSizeBiclique(string fpath)
{
    int maxSize = 0;

    ifstream file(fpath);
    string first;

    while (getline(file, first))
    {
        if (first.empty())
            continue;

        if (first.rfind("eps:", 0) == 0 || first.rfind("tau:", 0) == 0 || first.rfind("k", 0) == 0 || first.rfind("phi", 0) == 0)
            continue;

        string second;
        getline(file, second);

        int sizeBiclique = countSize(first) + countSize(second);
        maxSize = max(maxSize, sizeBiclique);
    }

    return maxSize;
}

double avgSizeBiclique(string fpath)
{
    int totSize = 0;
    int count = 0;

    ifstream file(fpath);
    string first;

    while (getline(file, first))
    {
        if (first.empty())
            continue;

        if (first.rfind("eps:", 0) == 0 || first.rfind("tau:", 0) == 0 || first.rfind("k", 0) == 0 || first.rfind("phi", 0) == 0)
            continue;

        string second;
        getline(file, second);

        int sizeBiclique = countSize(first) + countSize(second);
        totSize += sizeBiclique;
        count ++;
    }

    return static_cast<double>(totSize) / count;
}

void statsBicliques(string fpath)
{
    int maxSize = 0;
    int totSize = 0;
    int bicliqueCount = 0;

    ifstream file(fpath);
    string first;

    while (getline(file, first))
    {
        if (first.empty())
            continue;

        if (first.rfind("eps:", 0) == 0 || first.rfind("tau:", 0) == 0 || first.rfind("k", 0) == 0 || first.rfind("phi", 0) == 0)
            continue;

        string second;
        getline(file, second);

        int sizeBiclique = countSize(first) + countSize(second);
        maxSize = max(maxSize, sizeBiclique);
        totSize += sizeBiclique;

        bicliqueCount++;
    }

    double avgSize = static_cast<double>(totSize) / bicliqueCount;

    // print stats
    cout << fpath << endl;
    cout << "- Number of bicliques: " << bicliqueCount << endl;
    cout << "- Max size of bicliques: " << maxSize << endl;
    cout << "- Avg size of bicliques: " << avgSize << endl;
}

void statsBicliques(string ipath, ofstream &out)
{

    int maxSize = 0;
    int totSize = 0;
    int bicliqueCount = 0;

    ifstream file(ipath);
    string first;

    while (getline(file, first))
    {
        if (first.empty())
            continue;

        if (first.rfind("eps:", 0) == 0 || first.rfind("tau:", 0) == 0 || first.rfind("k", 0) == 0 || first.rfind("phi", 0) == 0)
            continue;

        string second;
        getline(file, second);

        int sizeBiclique = countSize(first) + countSize(second);
        maxSize = max(maxSize, sizeBiclique);
        totSize += sizeBiclique;

        bicliqueCount++;
    }

    double avgSize = 0;

    if (bicliqueCount > 0)
    {
        avgSize = static_cast<double>(totSize) / bicliqueCount;
    }

    // output stats
    out << bicliqueCount << " " << maxSize << " " << avgSize << endl;
}

void avgStats(string ipath, ofstream &out)
{
    ifstream file(ipath);

    int totBicliqueCount = 0;
    int totMaxSize = 0;
    double totAvgSize = 0;

    int n_run;
    int count;
    int maxSize;
    double avgSize;
    string line;
    getline(file, line);
    while (getline(file, line))
    {
        stringstream ss(line);
        ss >> n_run;
        ss >> count;
        ss >> maxSize;
        ss >> avgSize;

        totBicliqueCount += count;
        totMaxSize += maxSize;
        totAvgSize += avgSize;
    }

    double avgBicliqueCount = static_cast<double>(totBicliqueCount) / n_run;
    double avgMaxSize = static_cast<double>(totMaxSize) / n_run;
    double avgAvgSize = static_cast<double>(totAvgSize) / n_run;

    cout << "avgBicliqueCount: " << avgBicliqueCount << endl;
    cout << "avgMaxSize: " << avgMaxSize << endl;
    cout << "avgAvgSize: " << avgAvgSize << endl;

    out << avgBicliqueCount << " " << avgMaxSize << " " << avgAvgSize;
    
}

#endif