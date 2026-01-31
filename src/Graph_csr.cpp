#ifndef GRAPH_CSR_CPP
#define GRAPH_CSR_CPP

#include "include/Graph_csr.h"
#include "include/Utils.h"
#include "include/MinHashBall.h"

#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <algorithm>


Graph_csr::~Graph_csr()
{
    delete[] this->o_First;
    delete[] this->o_Target;
    delete[] this->o_degree;
    delete[] this->o_red_degree;

    delete[] this->counter; 

    if (this->directed)
    {
        delete[] this->i_First;
        delete[] this->i_Target;
        delete[] this->i_degree;
        delete[] this->o_red_degree;
    }

    delete[] this->balls;
    delete[] this->neigh;
}

// contructor

Graph_csr::Graph_csr(uint32_t N1, uint32_t N2, uint64_t M, bool isDirected, float Eps, int Tau)
{

    this->n1 = N1;
    this->n2 = N2;
    this->n = N1 + N2;
    this->m = M * (2 - (int)isDirected);
    this->directed = isDirected;

    this->o_First = new uint64_t[this->n];
    this->o_Target = new uint32_t[this->m];
    this->o_degree = new uint32_t[this->n];

    this->counter = new uint32_t[this->n1]; 

    if (isDirected)
    {
        this->i_First = new uint64_t[this->n];
        this->i_Target = new uint32_t[this->m];
        this->i_degree = new uint32_t[this->n];
    }

    this->del = new bool[this->n];
    this->eps = Eps;
    this->tau = Tau;
    this->neigh = new std::unordered_set<uint32_t>[this->n1];

    for (uint32_t i = 0; i < this->n; i++)
    {
        this->o_First[i] = 0;
        this->o_degree[i] = 0;

        if (i < n1)
        {
            this->counter[i] = 0;
        }

        if (isDirected)
        {
            this->i_First[i] = 0;
            this->i_degree[i] = 0;
        }

        this->del[i] = false;
    }

}

Graph_csr::Graph_csr(uint32_t N1, uint32_t N2, uint64_t M, bool isDirected, float Eps, int Tau, int K, float Phi, int n_hashes, Hash<uint32_t> **hash_functions)
    : Graph_csr(N1, N2, M, isDirected, Eps, Tau)
{
    
    this->o_red_degree = new uint32_t[this->n];

    if (isDirected)
    {
        this->i_red_degree = new uint32_t[this->n];
    }

    this->k = K;
    this->phi = Phi;

    this->balls = new MinHashBall[this->n];

    for (uint32_t i = 0; i < this->n; i++)
    {
        this->balls[i].insert(i);
        this->o_red_degree[i] = 0;

        if (isDirected)
        {
            this->i_red_degree[i] = 0;
        }
    }

    for (uint32_t u = 0; u < this->n; u++)
    {
        // this->balls[u] = MinHashBall(hash_functions, n_hashes);
        this->balls[u].init(hash_functions, n_hashes, u);
    }
}

Graph_csr *Graph_csr::from_file(std::string fname, bool isDirected, float eps, int tau)
{
    uint32_t n1;
    uint32_t n2;
    uint64_t m;
    uint32_t *edges = read_edges(fname, &n1, &n2, &m);
    Graph_csr *graph = Graph_csr::from_edges(edges, n1, n2, m, isDirected, eps, tau);
    delete[] edges;
    return graph;
}


Graph_csr *Graph_csr::from_edges(uint32_t *edges, uint32_t n1, uint32_t n2, uint64_t m, bool isDirected, float eps, int tau)
{
    Graph_csr *graph = new Graph_csr(n1, n2, m, isDirected, eps, tau);
    graph->process_edges(edges);
    return graph;
}

Graph_csr *Graph_csr::from_file(std::string fname, bool isDirected, float eps, int tau, int k, float phi, int n_hashes, Hash<uint32_t> **hash_functions)
{
    uint32_t n1;
    uint32_t n2;
    uint64_t m;
    uint32_t *edges = read_edges(fname, &n1, &n2, &m);
    Graph_csr *graph = Graph_csr::from_edges(edges, n1, n2, m, isDirected, eps, tau, k, phi, n_hashes, hash_functions);
    delete[] edges;
    return graph;
}


Graph_csr *Graph_csr::from_edges(uint32_t *edges, uint32_t n1, uint32_t n2, uint64_t m, bool isDirected, float eps, int tau, int k, float phi, int n_hashes, Hash<uint32_t> **hash_functions)
{
    Graph_csr *graph = new Graph_csr(n1, n2, m, isDirected, eps, tau, k, phi, n_hashes, hash_functions);
    graph->process_edges(edges);
    return graph;
}


void Graph_csr::process_edges(uint32_t *edges)
{
    uint32_t *max_indegree = new uint32_t[this->n];  // maximum in degree for each vertex; used to initialize i_First
    uint32_t *max_outdegree = new uint32_t[this->n]; // maximum out degree for each vertex; used to initialize o_First

    for (uint32_t i = 0; i < this->n; i++)
    {
        max_indegree[i] = 0;
        max_outdegree[i] = 0;
    }

    for (uint64_t i = 0; i < 2 * this->m / (2 - (int)this->directed); i += 2)
    {
        max_outdegree[edges[i]]++;
        if (this->directed)
            max_indegree[edges[i + 1]]++;
        else
            max_outdegree[edges[i + 1]]++;
    }

    // init i_First and o_First; no need to init the first vertex since it is already initialized
    for (uint32_t i = 1; i < this->n; i++)
    {
        this->o_First[i] = this->o_First[i - 1] + max_outdegree[i - 1];
        if (this->directed)
            this->i_First[i] = this->i_First[i - 1] + max_indegree[i - 1];
    }

    for (uint64_t i = 0; i < 2 * this->m / (2 - (int)this->directed); i += 2)
    {
        uint32_t u = edges[i];
        uint32_t v = edges[i + 1];
        this->insert_edge(u, v);
    }

    /*
    for (uint32_t i = 0; i < this->n; i++)
    {
        this->o_degree[i] = 0;
        if (this->directed)
            this->i_degree[i] = 0;
    }
    */

    delete[] max_indegree;
    delete[] max_outdegree;
}


void Graph_csr::insert_edge(uint32_t u, uint32_t v)
{
    // store (u, v) as out edge
    this->o_Target[this->o_First[u] + this->o_degree[u]] = v;
    this->o_degree[u]++;

    if (this->directed)
    {
        // store (u, v) as in edge
        this->i_Target[this->i_First[v] + this->i_degree[v]] = u;
        this->i_degree[v]++;
    }
    else
    {
        this->o_Target[this->o_First[v] + this->o_degree[v]] = u;
        this->o_degree[v]++;
    }
}


void Graph_csr::VReduce(bool apx)
{

    uint32_t *sim_deg = new uint32_t[this->n1];
    uint32_t *deg = new uint32_t[this->n];
    std::queue<uint32_t> Q;
    for (uint32_t u = 0; u < this->n; u++)
    {
        this->del[u] = false;
        deg[u] = this->o_degree[u];
        if (u < this->n1)
        {
            this->neigh[u] = this->SimNei(u, apx);
            sim_deg[u] = this->neigh[u].size();

            if (sim_deg[u] < this->tau - 1)
            {
                Q.push(u);
            }
        }

        if (deg[u] < this->tau)
            Q.push(u);
    }

    while (!Q.empty())
    {
        uint32_t u = Q.front();
        Q.pop();

        if (this->del[u]) continue;

        // ulteriore controllo
        if (deg[u] >= this->tau) continue;
        if (u < this->n1 && sim_deg[u] >= this->tau - 1) continue;

        this->del[u] = true;

        for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
        {
            uint32_t v = o_Target[i];
            if (this->del[v]) continue;

            if (--deg[v] < this->tau) Q.push(v);
        }

        if (u < this->n1){
            for (uint32_t v: this->neigh[u])
            {
                if (this->del[v]) continue;

                if (--sim_deg[v] < this->tau - 1) Q.push(v);
            }
        }
    }

    delete[] sim_deg;
    delete[] deg;
}


std::unordered_set<uint32_t> Graph_csr::SimNei(uint32_t u, bool apx)
{
    if (!apx)
    {
        // jaccard 1 hop
        std::unordered_set<uint32_t> sim_nei;
        
        if (this->del[u]) return sim_nei;

        std::vector<uint32_t> touched;

        for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
        {
            uint32_t v = o_Target[i];

            for (uint32_t j = o_First[v]; j < o_First[v] + o_degree[v]; j++)
            {
                uint32_t w = o_Target[j];

                if (u == w || this->del[w]) continue;

                if (this->counter[w] == 0) touched.push_back(w);
                
                this->counter[w]++;
            }
        }

        for (uint32_t v : touched)
        {
            double sim = 
                static_cast<double>(this->counter[v]) / 
                (o_degree[u] + o_degree[v] - this->counter[v]);
            
                if (sim >= this->eps) sim_nei.insert(v);
            
                this->counter[v] = 0;
        }

        return sim_nei;
    }
    else
    {
        //jaccard 2 hop
        std::unordered_set<uint32_t> sim_nei;

        for (uint32_t v = 0; v < this->n1; v++)
        {
            if (MinHashBall::similarity((&this->balls[u]), (&this->balls[v])) >= this->eps) sim_nei.insert(v);
        }

        return sim_nei;
    }
}

void Graph_csr::MSBE(std::string fname, bool apx)
{

    std::string fpath = "./dataset/results/bicliques/" + fname + (apx? "_2_hop": "_1_hop") + ".txt";

    std::ofstream out(fpath, std::ios::trunc);

    out << "eps: " << this->eps << endl;
    out << "tau: " << this->tau << endl;

    if (apx)
    {
        out << "k: " << this->k << endl;
        out << "phi: " << this->phi << endl;
    }

    for (uint32_t u = 0; u < this->n1; u++)
    {
        if (this->del[u]) continue;

        std::unordered_set<uint32_t> C_l = {u};
        std::unordered_set<uint32_t> C_r;
        std::unordered_set<uint32_t> P_l;
        std::unordered_set<uint32_t> Q_l;

        for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
        {
            uint32_t v = o_Target[i];
            if (this->del[v]) continue;
            C_r.insert(v);
        }

        for (uint32_t v: this->neigh[u])
        {
            if (v > u) P_l.insert(v);
            else if (v < u) Q_l.insert(v);
        }

        this->Enum(C_l, C_r, P_l, Q_l, fpath);
    }
}

std::unordered_set<uint32_t> set_union(const std::unordered_set<uint32_t> &a, const std::unordered_set<uint32_t> &b)
{
    std::unordered_set<uint32_t> union_set = a;
    union_set.insert(b.begin(), b.end());
    return union_set;
}

std::unordered_set<uint32_t> set_intersection(const std::unordered_set<uint32_t> &a, const std::unordered_set<uint32_t> &b)
{
    const auto& small = (a.size() < b.size()) ? a : b;
    const auto& large = (a.size() < b.size()) ? b : a;

    std::unordered_set<uint32_t> result;
    result.reserve(small.size());

    for (const auto& x : small) {
        if (large.find(x) != large.end()) {
            result.insert(x);
        }
    }
    return result;
}

void outputBiclique(std::unordered_set<uint32_t> &C_l, std::unordered_set<uint32_t> &C_r, std::string fpath)
{
    std::ofstream out(fpath, std::ios::app);
    
    std::vector<uint32_t> L(C_l.begin(), C_l.end());
    std::vector<uint32_t> R(C_r.begin(), C_r.end());
    std::sort(L.begin(), L.end());
    std::sort(R.begin(), R.end());

    out << "\n";

    for(auto& x: L)
    {
        out << x << " ";
    }
    out << "\n";

    for(auto& x: R)
    {
        out << x << " ";
    }
    out << "\n";

}

void Graph_csr::Enum(std::unordered_set<uint32_t> C_l, std::unordered_set<uint32_t> C_r, std::unordered_set<uint32_t> P_l, std::unordered_set<uint32_t> Q_l, std::string fpath)
{

    bool isMaximal = true;
    for (uint32_t u: P_l)
    {
        // check whether u structurally connects to all nodes in C_r
        bool extends = true;
        for (uint32_t v: C_r)
        {
            bool found = false;
            for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
            {
                if (o_Target[i] == v)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                extends = false;
                break;
            }
        }

        // if u extends C_r, biclique is NOT maximal
        if (extends) 
        {
            isMaximal = false;
        }
    }

    for (uint32_t u: Q_l)
    {
        // check whether u structurally connects to all nodes in C_r
        bool extends = true;
        for (uint32_t v: C_r)
        {
            bool found = false;
            for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
            {
                if (o_Target[i] == v)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                extends = false;
                break;
            }
        }

        // if u extends C_r, biclique is NOT maximal
        if (extends) 
        {
            isMaximal = false;
            
            //early termination
            bool sim2PL = true;
            for (auto v: P_l)
            {
                if (this->neigh[u].find(u) == this->neigh[u].end())
                {
                    sim2PL = false;
                }
            }
            if (sim2PL) return;
        }
    }

    if (isMaximal && C_r.size() >= this->tau && C_l.size() >= this->tau)
    {

        outputBiclique(C_l, C_r, fpath);
    }

   
    std::vector P_snapshot(P_l.begin(), P_l.end());
    for (uint32_t u: P_snapshot)
    {
        
        P_l.erase(u);

        std::unordered_set<uint32_t> C_l1 = C_l;
        std::unordered_set<uint32_t> C_r1;
        std::unordered_set<uint32_t> P_l1;
        std::unordered_set<uint32_t> Q_l1;
        C_l1.insert(u);
        std::unordered_set<uint32_t> N_u; 
        for (uint32_t i = o_First[u]; i < o_First[u] + o_degree[u]; i++)
        {
            N_u.insert(o_Target[i]);
        }
        C_r1 = set_intersection(C_r, N_u);
        P_l1 = set_intersection(P_l, this->neigh[u]);
        Q_l1 = set_intersection(Q_l, this->neigh[u]);

        if (C_l1.size() + P_l1.size() >= this->tau && C_r1.size() >= this->tau)
        {
            Enum(C_l1, C_r1, P_l1, Q_l1, fpath);
        }

        Q_l.insert(u);

    }
    
}

int Graph_csr::update(uint32_t u, uint32_t v)
{
    int n_merge = 0;

    // increment the red degree of u
    this->o_red_degree[u]++;

    // increment the degree of u
    this->o_degree[u]++;

    // insert v into the ball of u
    this->balls[u].insert(v);
    n_merge++;

    // push the ball of radius 1 of v to the ball of radius 2 of u
    this->balls[u].push(&this->balls[v]);
    n_merge++;

    // check if the red degree of u is greater than the threshold
    float threshold = this->phi * (this->o_degree[u] - this->o_red_degree[u]);
    if (this->o_red_degree[u] >= threshold)
    {
        // reset the red degree of u
        this->o_red_degree[u] = 0;
        // propagate the verites of u's ball of radius 1 to the ball of radius 2 of all its neighbours
        n_merge += this->propagate(u);
    }
    else
    {
        n_merge += min((uint32_t)this->k, !this->directed ? this->o_degree[u] : this->i_degree[u]);
        for (int i = 0; i < this->k; i++)
        {
            if (!this->directed)
            {
                uint64_t rand_index = this->o_First[u] + (rand() % this->o_degree[u]);
                uint32_t x = this->o_Target[rand_index];
                this->balls[x].push(&this->balls[u]);
            }
            else
            {
                if (this->i_degree[u] == 0)
                    break;
                uint64_t rand_index = this->i_First[u] + (rand() % this->i_degree[u]);
                uint32_t x = this->i_Target[rand_index];
                this->balls[x].push(&this->balls[u]);
            }
        }
    }

    if (!this->directed)
    {
        this->o_red_degree[v]++;
        this->o_degree[v]++;
        this->balls[v].insert(u);
        n_merge++;
        this->balls[v].push(&this->balls[u]);
        n_merge++;

        threshold = this->phi * (this->o_degree[v] - this->o_red_degree[v]);
        if (this->o_red_degree[v] >= threshold)
        {
            this->o_red_degree[v] = 0;
            n_merge += this->propagate(v);
        }
        else
        {
            n_merge += min((uint32_t)this->k, this->o_degree[v]);
            for (int i = 0; i < k; i++)
            {
                uint64_t rand_index = this->o_First[v] + (rand() % this->o_degree[v]);
                uint32_t x = o_Target[rand_index];
                this->balls[x].push(&this->balls[v]);
            }
        }
    }
    else
    {
        this->i_degree[v]++;
    }

    return n_merge;
}

/**
 * This method is used to propagate the vertices of the ball of radius 1 of the given vertex u to the ball of radius 2 of all its neighbours.
 * @param u: vertex whose ball of radius 1 is to be propagated to the ball of radius 2 of all its neighbours.
 * @return: void
 */
int Graph_csr::propagate(uint32_t u)
{
    if (this->directed)
    {
        for (uint64_t i = this->i_First[u]; i < this->i_First[u] + this->i_degree[u]; i++)
        {
            uint32_t v = this->i_Target[i];
            this->balls[v].push(&this->balls[u]);
        }
        return this->i_degree[u];
    }
    else
    {
        for (uint64_t i = this->o_First[u]; i < this->o_First[u] + this->o_degree[u]; i++)
        {
            uint32_t v = this->o_Target[i];
            this->balls[v].push(&this->balls[u]);
        }
        return this->o_degree[u];
    }
}

void Graph_csr::flush_graph()
{
    for (uint32_t i = 0; i < this->n; i++)
    {
        this->o_degree[i] = 0;
        if (this->directed)
            this->i_degree[i] = 0;

        this->o_red_degree[i] = 0;
        this->balls[i].flush(i);
    }
}

void Graph_csr::printDel()
{
    uint32_t pruned = 0;
    uint32_t remaining = 0;
    for (uint32_t u = 0; u < this->n; u++)
    {
        if (this->del[u])
        {
            pruned++;
        }else{
            remaining++;
        }
    }
    std::cout << "- pruned nodes: " << pruned << "\n- remaining nodes: " << remaining << std::endl;
}

#endif