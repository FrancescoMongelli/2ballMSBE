#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <string>
#include <fstream>

int numberOfBicliques(std::string fpath);
int countSize(std::string line);
int maxSizeBiclique(std::string fpath);
double avgSizeBiclique(std::string fpath);
void statsBicliques(std::string fpath);
void statsBicliques(std::string ipath, std::ofstream &out);
void avgStats(std::string ipath, std::ofstream &out);

#endif

