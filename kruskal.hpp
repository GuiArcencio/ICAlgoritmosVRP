#ifndef KRUSKAL_HPP
#define KRUSKAL_HPP

#include <vector>

#include "unionfind.hpp"

namespace kruskal
{
    struct edge
    {
        int u, v;
        double length;
        bool starting;
    };

    std::vector<edge> generateMSG(int n, int K, std::vector<edge>& edges, double* cost, double* max_length);
}

#endif