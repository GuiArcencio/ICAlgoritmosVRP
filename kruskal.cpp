#include "kruskal.hpp"

#include <algorithm>
#include <stdio.h>

std::vector<kruskal::edge> kruskal::generateMSG(int n, int K, std::vector<kruskal::edge>& edges, double* cost, double* max_length)
{
    int extra_edges = 0;
    *cost = 0.0;
    *max_length = 0.0;
    UnionFind components(n);

    std::vector<kruskal::edge> msg;
    msg.reserve(n - 1 + K);

    // Ordena colocando as arestas iniciais no começo do vetor
    std::sort(edges.begin(), edges.end(), [] (const kruskal::edge& a, const kruskal::edge& b) -> bool {
        if (a.starting && !b.starting) 
            return true;
        else if (!a.starting && b.starting)
            return false;
        else return a.length < b.length;
    });

    for (auto e : edges)
    {
        // Adding preselected edges
        if (e.starting)
        {
            msg.push_back(e);
            *cost += e.length;
            if (components.find(e.u) == components.find(e.v))
                extra_edges++;
            else
                components.unite(e.u, e.v);

            if (e.length > *max_length)
                *max_length = e.length;
        }
        else
        {
            if (components.find(e.u) != components.find(e.v))
            {
                msg.push_back(e);
                *cost += e.length;
                components.unite(e.u, e.v);

                if (e.length > *max_length)
                    *max_length = e.length;
            }
            else if (extra_edges < K)
            {
                msg.push_back(e);
                *cost += e.length;
                extra_edges++;

                if (e.length > *max_length)
                    *max_length = e.length;
            }
        }
    }

    return msg;
}