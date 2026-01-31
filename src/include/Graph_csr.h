#ifndef GRAPH_CSR_H
#define GRAPH_CSR_H

#include "../Hash.cpp"
#include "./MinHashBall.h"

#include <cstdint>
#include <vector>
#include <unordered_set>
#include <string>

class Graph_csr
{

private:
    /// Number of vertices in the graph
    uint32_t n1;
    uint32_t n2;
    uint32_t n;

    /// Number of edges in the graph
    uint64_t m;

    /// Boolean value indicating whether the graph is directed or not
    bool directed;

    uint64_t *o_First;          // index, in o_Target, where the outedges of vertex v starts. Size n
    uint32_t *o_Target;         // second enpoint for each edge. Size m
    uint32_t *o_degree;         // current out-degree for each vertex v. Size n
    uint32_t *o_red_degree;

    uint32_t *counter;

    uint64_t *i_First;          // index, in i_Target, where the inedges of vertex v starts. Size n
    uint32_t *i_Target;         // second enpoint for each edge. Size m
    uint32_t *i_degree;         // current in-degree for each vertex v. Size n
    uint32_t *i_red_degree;

    bool *del;                  // boolean value indicating if the vertex v is pruned. Size n
    float eps;                 // precision requirement for the sim function between members of one community
    int tau;               // minimum size of the communities
    int k;
    float phi;
    MinHashBall *balls;          // 2 balls of the vertices. Size n
    std::unordered_set<uint32_t> *neigh;    // similar neighborhood. Size n

    void process_edges(uint32_t *edges);

public:

    /**
     * This is the default destructor of the Graph_csr class.
     */
    ~Graph_csr();

    /**
     * This is the default constructor of the Graph_csr class.
     * @param N: number of vertices in the graph.
     * @param M: number of edges in the graph.
     * @param isDirected: boolean value indicating whether the graph is directed or not.
     * @param Eps: precision requirement for the sim function between members of one community
     * @param Tau: minimum size of the communities
     */
    Graph_csr(uint32_t N1, uint32_t N2, uint64_t M, bool isDirected, float Eps, int Tau);

    Graph_csr(uint32_t N1, uint32_t N2, uint64_t M, bool isDirected, float Eps, int Tau, int K, float Phi, int n_hashes, Hash<uint32_t> **hash_functions);

    static Graph_csr *from_file(std::string fname, bool isDirected, float eps, int tau);

    static Graph_csr *from_edges(uint32_t *edges, uint32_t n1, uint32_t n2, uint64_t m, bool isDirected, float eps, int tau);

    static Graph_csr *from_file(std::string fname, bool isDirected, float eps, int tau, int k, float phi, int n_hashes, Hash<uint32_t> **hash_functions);

    static Graph_csr *from_edges(uint32_t *edges, uint32_t n1, uint32_t n2, uint64_t m, bool isDirected, float eps, int tau, int k, float phi, int n_hashes, Hash<uint32_t> **hash_functions);

    /**
     * This method inserts the edge (u, v) into the graph.
     * @param u: first endpoint of the edge to insert.
     * @param v: second endpoint of the edge to insert.
     * @return: void
     */
    void insert_edge(uint32_t u, uint32_t v);

    /**
     * This method popolates the function del, which exludes the vertices which can't create the cliques.
     * @param 
     */
    void VReduce(bool apx = false);

    /**
     * This method returns the neighbours of the vertex u that are similar to u
     * @param u: vertex
     * @return: the minilar neighborhood of the vertex u
     */
    std::unordered_set<uint32_t> SimNei(uint32_t u, bool apx = false);

    /**
     * Core of the alghorithm of MSBE
     * prints to a file the communities in the format:
     * leftnode1 leftnode2...
     * rightnode1 rightnode2...
     */
    void MSBE(std::string fname, bool apx = false);

    void Enum(std::unordered_set<uint32_t> C_l, std::unordered_set<uint32_t> C_r, std::unordered_set<uint32_t> P_l, std::unordered_set<uint32_t> Q_l, std::string fpath);

    int update(uint32_t u, uint32_t v);

    int propagate(uint32_t u);

    void flush_graph();

    // DA ELIMINAREEEEEEEEEEEEEEEE
    void printDel();

};

#endif