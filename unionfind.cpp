#include "unionfind.hpp"

UnionFind::UnionFind(unsigned int n)
{
    this->parent = new unsigned int[n];
    this->rank = new unsigned int[n];

    for (unsigned int i = 0; i < n; i++)
    {
        this->parent[i] = i;
        this->rank[i] = 0;
    }
}

UnionFind::~UnionFind()
{
    delete[] this->parent;
    delete[] this->rank;
}

unsigned int UnionFind::find(unsigned int x)
{
    if (x != parent[x])
        parent[x] = find(parent[x]);
    return parent[x];
}

void UnionFind::unite(unsigned int x, unsigned int y)
{
    unsigned int rx = find(x);
    unsigned int ry = find(y);

    if (rx == ry) return;
    if (rank[rx] > rank[ry])
        parent[ry] = rx;
    else
    {
        parent[rx] = ry;
        if (rank[rx] == rank[ry]) 
            rank[ry] += 1;
    }
}