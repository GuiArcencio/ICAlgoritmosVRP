#ifndef SPANNING_COVER_HPP
#define SPANNING_COVER_HPP

#include <vector>
#include <list>

#include "gurobi_c++.h"

namespace spi
{

struct edge
{
    int u, v;
    double d;
};

double getDistance(GRBVar** x, int u, int v);

double getDelta(GRBVar** vars, const std::vector<std::list<edge>>& graph, edge starting, edge* e_prime);

double getLowerBound(GRBVar** vars, const std::vector<edge>& ordered, const std::list<edge>& starting, int n, int k);

}

#endif