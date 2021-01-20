#ifndef UNIONFIND_HPP
#define UNIONFIND_HPP

class UnionFind
{
    private:
        unsigned int *parent;
        unsigned int *rank;

    public:
        UnionFind(unsigned int n);
        ~UnionFind();

        unsigned int find(unsigned int x);
        void unite(unsigned int x, unsigned int y);
};

#endif