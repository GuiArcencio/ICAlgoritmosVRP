#ifndef KARGER_HPP
#define KARGER_HPP

#include <list>
#include <random>
#include <vector>
#include <map>

#include "unionfind.hpp"

namespace karger
{

struct edge
{
    unsigned int u;
    unsigned int v;
    double weight;
    double accum_weight;
    bool chosen;
};

class EdgeVector
{
    private:
        std::vector<edge> start_edges;
        unsigned int num_v;
        std::default_random_engine generator;

    public:
        EdgeVector(unsigned int num_vertices, unsigned int seed);
        ~EdgeVector();

        void add_edge(unsigned int u, unsigned int v, double weight);
        std::vector<std::list<int>> randomCut(const double* demands, std::vector<double>& sum_of_demands, int K);
        void clear_edges();
};

}

#endif