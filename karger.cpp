#include "karger.hpp"

karger::EdgeVector::EdgeVector(unsigned int num_vertices, unsigned int seed):
    num_v(num_vertices), generator(seed) 
{
    start_edges.reserve(num_v * num_v);
}
karger::EdgeVector::~EdgeVector() {}

void karger::EdgeVector::add_edge(unsigned int u, unsigned int v, double weight)
{
    if (weight > 0.0001)
    {
        edge e;
        e.u = u;
        e.v = v;
        e.chosen = false;
        e.weight = weight;
        if (start_edges.empty())
            e.accum_weight = weight;
        else
            e.accum_weight = weight + start_edges[start_edges.size() - 1].accum_weight;
            
        start_edges.push_back(e);
    }
}

std::vector<std::list<int>> karger::EdgeVector::randomCut(const double* demands, std::vector<double>& sum_of_demands, int K)
{
    UnionFind merges(num_v);
    int current_n = num_v - 1; // excludes depot
    double max_weight = start_edges[start_edges.size() - 1].accum_weight;
    int remaining_edges = start_edges.size();

    while (current_n > K && remaining_edges > 0)
    {
        unsigned int uv;

        std::uniform_real_distribution<double> prob_dist(0.0, max_weight);
        double total_value = prob_dist(generator);

        unsigned int left = 0;
        unsigned int right = start_edges.size() - 1;
        bool found = false;

        while (!found)
        {
            unsigned int middle = left + (right - left)/2;

            // Found
            if ( total_value < start_edges[middle].accum_weight && 
                total_value >= (start_edges[middle].accum_weight - start_edges[middle].weight) )
            {
                uv = middle;
                found = true;
            }
            else if (total_value >= start_edges[middle].accum_weight)
                left = middle + 1;
            else
                right = middle - 1;
        }

        // uv is now the index to the chosen edge

        if (!start_edges[uv].chosen)
        {
            if (merges.find(start_edges[uv].u) != merges.find(start_edges[uv].v))
            {
                merges.unite(start_edges[uv].u, start_edges[uv].v);
                --current_n;
            }
            start_edges[uv].chosen = true;
            remaining_edges--;
        }
    }

    sum_of_demands.resize(K, 0.0);
    std::vector<std::list<int>> cuts(K);

    std::map<int, int> vertex_to_index;
    int next_index = 0;
    for (int i = 1; i < num_v; i++)
        if (merges.find(i) == i)
        {
            vertex_to_index[i] = next_index;
            ++next_index;
        }

    for (int i = 1; i < num_v; i++)
    {
        int index = vertex_to_index[merges.find(i)];
        if (index < K)
        {
            cuts[index].push_back(i);
            sum_of_demands[index] += demands[i];
        }
    }

    for (edge& e : start_edges)
        e.chosen = false;

    return cuts;
}

void karger::EdgeVector::clear_edges()
{
    start_edges.clear();
    start_edges.reserve(num_v * num_v);
}