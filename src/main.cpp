#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "include/Graph_csr.h"
#include "include/Utils.h"
#include "include/Experiments.h"
#include "./Hash.cpp"

using namespace std;

void classicMSBE(string fname, bool isDirected, float eps, int tau)
{
    cerr << fname << endl
         << flush;

    Graph_csr *g = Graph_csr::from_file(fname, isDirected, eps, tau);

    filesystem::path p(fname);
    string dataname = p.stem().string();
    
    g->VReduce(false);
    cerr << "VReduce done\n" << flush;
    g->printDel();
    g->MSBE(dataname, false);
}

void apxMSBE(string fname, bool isDirected, float eps, int tau, int k, float phi, int n_hashes)
{
    cerr << fname << endl
         << flush;

    TabulationHash<uint32_t> **hashes = new TabulationHash<uint32_t> *[n_hashes]; 
    for (int i = 0; i < n_hashes; i++)
        hashes[i] = new TabulationHash<uint32_t>();

    Graph_csr *g = Graph_csr::from_file(fname, isDirected, eps, tau, k, phi, n_hashes, (Hash<uint32_t> **)hashes);
    cerr << "Graph reading done\n" << flush;

    uint32_t n1;
    uint32_t n2;
    uint64_t m;
    uint32_t *edges = read_edges(fname, &n1, &n2, &m);

    g->flush_graph();
    for(int i=0; i< 2*m; i+=2){
        g->update(edges[i], edges[i+1]);
    }
    cerr << "2 balls updated\n" << flush;

    filesystem::path p(fname);
    string dataname = p.stem().string();
    
    g->VReduce(true);
    cerr << "VReduce done\n" << flush;
    g->printDel();
    g->MSBE(dataname, true);

    delete[] edges;
    delete g;

    for (int i = 0; i < n_hashes; i++)
        delete hashes[i];
    delete[] hashes;
}

void compareBicliques(string dataname, ofstream &out, int K = 48, int r = 6) //K=128 K=8,r=4
{
    double D = 1.0;

    vector<Biclique> bicliques;

    auto A = parseBicliqueFile("./dataset/results/bicliques/" + dataname + "_1_hop.txt", false);
    auto B = parseBicliqueFile("./dataset/results/bicliques/" + dataname + "_2_hop.txt", true);

    int nA = A.size();
    int nB = B.size();

    bicliques.insert(bicliques.end(), A.begin(), A.end());
    bicliques.insert(bicliques.end(), B.begin(), B.end());


    if (nA == 0 || nB == 0) {
        cerr << "Warning: one file has zero bicliques\n";
        out << D << endl;
        return;
    }

    vector<Hash<uint32_t>*> hashes;
    hashes.reserve(K);

    for (int i = 0; i < K; i++)
        hashes.push_back(new TabulationHash<uint32_t>());

    int b = K / r;

    auto candidates = runLSH(bicliques, hashes, r, b);


    D = dissimilarityIndex(candidates, nA, nB);

    for (auto h : hashes)
        delete h;

    cout << "Dissimilarity: " << D*100 << "%\n";
    out << D << endl;
}

int main(int argc, const char *argv[])
{

    string usage = "./main <dataset> <isDirected> <epsilon> <tau> [<k> <phi> <n_hashes> [<n_run>]]";
    if (argc < 5)
    {
        cout << usage << endl;
        return 0;
    }

    string filename = argv[1];
    bool isDirected = (bool)atoi(argv[2]);
    float eps = atof(argv[3]);
    int tau = atoi(argv[4]);

    if (argc == 5)
    {
        classicMSBE(filename, isDirected, eps, tau);

        filesystem::path p(filename);
        string dataname = p.stem().string();
        string ipath = "./dataset/results/bicliques/" + dataname + "_1_hop.txt";
        string opath = "./dataset/results/stats/stats_" + dataname + "_1_hop.txt";

        ofstream out(opath, ios::trunc);
        out << "n_bicliques max_size avg_size\n";

        statsBicliques(ipath, out);
    }
    else if (argc == 8)
    {
        int k = atoi(argv[5]);
        float phi = atof(argv[6]);
        int n_hashes = atoi(argv[7]);

        apxMSBE(filename, isDirected, eps, tau, k, phi, n_hashes);
    }
    else if (argc == 9)
    {
        
        int k = atoi(argv[5]);
        float phi = atof(argv[6]);
        int n_hashes = atoi(argv[7]);
        int n_runs = atoi(argv[8]);

        
        filesystem::path p(filename);
        string dataname = p.stem().string();
        string ipath = "./dataset/results/bicliques/" + dataname + "_2_hop.txt";
        string opath = "./dataset/results/stats/stats_" + dataname + "_2_hop.txt";
        string dispath = "./dataset/results/dissimilarity/dis_" + dataname + ".txt";

        ofstream out(opath, ios::trunc);
        out << "n_run n_bicliques max_size avg_size\n";

        ofstream disout(dispath, ios::trunc);
        disout << "n_run dissimilarity\n";

        for (int i = 0; i < n_runs; i++)
        {
            apxMSBE(filename, isDirected, eps, tau, k, phi, n_hashes);

            out << i+1 << " ";
            disout << i+1 << " ";

            statsBicliques(ipath, out);

            compareBicliques(dataname, disout);

            cout << "run " << i+1 << " completed\n";
        }
        out.close();

        ofstream statsOut("./dataset/results/stats/avg_stats_" + dataname + "_2_hop.txt", ios::trunc);
        statsOut << "n_bicliques max_size avg_size\n";
        avgStats(opath, statsOut);

        ofstream avgdisout("./dataset/results/dissimilarity/avg_dis_" + dataname + ".txt", ios::trunc);
        avgdisout << "n_run avg_dissimilarity std_deviation\n";
        avgDissimilarity(dispath, avgdisout);
    }
    else
    {
        cout << usage << endl;
    }

}