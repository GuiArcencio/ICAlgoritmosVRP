#include "spanningcover.hpp"

#include <stack>
#include <cstdio>

#include "unionfind.hpp"

double spi::getDistance(GRBVar** x, int u, int v) 
{
    if (u > v)
        return x[u][v].get(GRB_DoubleAttr_Obj);
    else if (v > u)
        return x[v][u].get(GRB_DoubleAttr_Obj);

    return -1;
}

double spi::getDelta(GRBVar** vars, const std::vector<std::list<edge>>& graph, edge starting, edge* e_prime)
{
    int n = graph.size();

    // Running a DFS to find the way back from e.v to the warehouse
    // (We suppose starting.u is the warehouse and starting.v is the next vertex)
    bool found_way = false;
    bool* discovered = new bool[n];
    int* parent = new int[n];
    for (int i = 0; i < n; i++) 
    {
        discovered[i] = false;
        parent[i] = -1;
    }

    std::stack<int> S;
    S.push(starting.v);

    while (!found_way)
    {
        int current = S.top();
        S.pop();

        if (current == 0)
            found_way = true;
        else if (!discovered[current])
        {
            discovered[current] = true;
            for (auto e : graph[current])
            {
                S.push(e.v);
                if(parent[e.v] == -1) parent[e.v] = current;
            }
        }
    }

    // Now looking for the most expensive edge on the path
    int current = parent[0];
    int next = parent[current];
    double max_cost = spi::getDistance(vars, current, next);
    *e_prime = { current, next, max_cost };

    current = next;
    next = parent[next];
    while (current != starting.v)
    {
        if (spi::getDistance(vars, current, next) > max_cost)
        {
            max_cost = spi::getDistance(vars, current, next);
            *e_prime = { current, next, max_cost };
        }

        current = next;
        next = parent[next];
    }

    delete[] discovered;
    delete[] parent;

    return starting.d - max_cost;
}

double spi::getLowerBound(GRBVar** vars, const std::vector<spi::edge>& ordered, const std::list<spi::edge>& starting, int n, int k)
{
    UnionFind connectivity(n);
    std::vector<std::list<spi::edge>> spanning_subgraph(n);
    int warehouse_edges = 0;
    int total_edges = 0;
    double total_cost = 0.0;

    // First, build a spanning subgraph which is minimal among those with the starting edges
    for (auto e : starting)
    {
        connectivity.unite(e.u, e.v);
        total_cost += e.d;
        total_edges++;

        spanning_subgraph[e.u].push_back({ e.u, e.v, e.d });
        spanning_subgraph[e.v].push_back({ e.v, e.u, e.d });

        if (e.u == 0 || e.v == 0)
            warehouse_edges++;
    }

    for (auto e : ordered)
    {
        if (connectivity.find(e.u) != connectivity.find(e.v))
        {
            connectivity.unite(e.u, e.v);
            total_cost += e.d;
            total_edges++;

            spanning_subgraph[e.u].push_back({ e.u, e.v, e.d });
            spanning_subgraph[e.v].push_back({ e.v, e.u, e.d });
            
            if (e.u == 0 || e.v == 0)
                warehouse_edges++;
        }
    }
    
    // Making exchanges until we have 2*k warehouse edges
    while (warehouse_edges < 2*k)
    {
        // Looking for the minimum delta
        edge to_be_added, to_be_removed;
        double min_delta = -1;

        for (int v = 1; v < n; v++)
        {
            bool already_in_ssg = false;
            for (auto it = spanning_subgraph[0].begin(); it != spanning_subgraph[0].end() && !already_in_ssg; ++it)
                if (it->v == v) already_in_ssg = true;

            if (!already_in_ssg)
            {
                edge e_prime;
                double delta = getDelta(vars, spanning_subgraph, { 0, v, getDistance(vars, 0, v) }, &e_prime);

                if (min_delta == -1 || delta < min_delta)
                {
                    min_delta = delta;
                    to_be_added = { 0, v, getDistance(vars, 0, v) };
                    to_be_removed = e_prime;
                }
            }
        }

        spanning_subgraph[to_be_removed.u].remove_if([to_be_removed] (const edge& e) {
            return e.v == to_be_removed.v;
        });
        spanning_subgraph[to_be_removed.v].remove_if([to_be_removed] (const edge& e) {
            return e.v == to_be_removed.u;
        });
        spanning_subgraph[to_be_added.u].push_back({ to_be_added.u, to_be_added.v, to_be_added.d });
        spanning_subgraph[to_be_added.v].push_back({ to_be_added.v, to_be_added.u, to_be_added.d });

        total_cost += min_delta;
        warehouse_edges++;
    }

    // Adding redundant edges up to n-1+k total
    for (int i = 0; i < ordered.size() && total_edges < n-1+k; i++)
    {
        bool already_in_ssg = false;
        for (auto it = spanning_subgraph[ordered[i].u].begin(); it != spanning_subgraph[ordered[i].u].end() && !already_in_ssg; ++it)
            if (it->v == ordered[i].v) already_in_ssg = true;

        if (!already_in_ssg)
        {
            total_cost += ordered[i].d;
            total_edges++;
        }
    }

    return total_cost;
}