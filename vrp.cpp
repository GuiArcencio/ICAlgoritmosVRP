#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <list>

#include "gurobi_c++.h"
#include "karger.hpp"

struct Point 
{
    double x;
    double y;
    double d;
};

double dist(const Point& a, const Point& b);

std::vector<Point> getPointsFromFile(std::string filename, int *V, double *C);
double** getHeuristicSol(std::string filename, int N, int V);
void writeSolution(GRBVar** x, int N, int V, const double obj, const int opt, const std::string& filename);

class subtourelim: public GRBCallback
{
    public:
        GRBVar** x;
        int N, V;
        double C;
        double* demands;
        karger::EdgeVector cut_generator;
        subtourelim(GRBVar** x, double* demands, int N, int V, double C):
            x(x), N(N), V(V), C(C), demands(demands), cut_generator(N, 1) {};

    protected:
        void callback()
        {
            if (where == GRB_CB_MIPSOL)
            {
                std::vector<bool> seen;
                for (int i = 0; i < N; i++)
                    seen.push_back(false);
                int start, current = 0;
                double route_demand;
                bool found_way;
                
                // Sweeping through depot routes
                do {
                    start = 0;
                    for (int i = 1; i < N; i++)
                        if (!seen[i] && getSolution(x[i][0]) > 0.5)
                        {
                            start = i;
                            break;
                        }

                    if (start != 0)
                    {
                        std::list<int> tour;
                        tour.push_back(start);
                        route_demand = demands[start];
                        found_way = false;
                        current = start;
                        seen[start] = true;

                        do {
                            found_way = false;
                            for (int j = 1; j < N; j++)
                                {
                                    if (current != j && !seen[j] && getSolution(x[std::max(current, j)][std::min(current, j)]) > 0.5)
                                    {
                                        found_way = true;
                                        tour.push_back(j);
                                        current = j;
                                        seen[j] = true;
                                        route_demand += demands[j];
                                        break;
                                    }
                                }
                        } while(found_way);

                        // if the route exceeds capacity
                        if (route_demand > C)
                        {
                            GRBLinExpr c = 0.0;
                            double ub = std::ceil(route_demand / C);
                            for (auto i = tour.begin(); i != tour.end(); ++i)
                                for (auto j = std::next(i); j != tour.end(); ++j)
                                    c += x[std::max(*i,*j)][std::min(*i,*j)];

                            addLazy(c, GRB_LESS_EQUAL, tour.size() - ub);
                        }
                    }
                    
                } while (start != 0);

                // Sweeping through non-depot routes
                for (start = 1; start < N; start++)
                {
                    if (!seen[start])
                    {
                        std::list<int> tour;
                        tour.push_back(start);
                        route_demand = demands[start];
                        seen[start] = true;
                        current = start;

                        do {
                            found_way = false;
                            for (int j = 1; j < N; j++)
                                {
                                    if (current != j && !seen[j] && getSolution(x[std::max(current, j)][std::min(current, j)]) > 0.5)
                                    {
                                        found_way = true;
                                        tour.push_back(j);
                                        current = j;
                                        seen[j] = true;
                                        route_demand += demands[j];
                                        break;
                                    }
                                }
                        } while (found_way);

                        GRBLinExpr c = 0.0;
                        double ub = std::ceil(route_demand / C);
                        for (auto i = tour.begin(); i != tour.end(); ++i)
                            for (auto j = std::next(i); j != tour.end(); ++j)
                                c += x[std::max(*i,*j)][std::min(*i,*j)];

                        addLazy(c, GRB_LESS_EQUAL, tour.size() - ub);
                    }
                } 

            }
            else if (where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
            {
                std::vector<double> sum_of_demands;

                for (int i = 1; i < N; i++)
                    for (int j = 1; j < i; j++)
                        cut_generator.add_edge(i, j, getNodeRel(x[i][j]));

                // Run Karger 
                for (int count = 0; count < 10 * N; count++)
                {
                    sum_of_demands.clear();
                    auto cuts = cut_generator.randomCut(demands, sum_of_demands, V);

                    // Going through each cut
                    for (int cut_i = 0; cut_i < cuts.size(); ++cut_i)
                    {
                        double sum_demands = sum_of_demands[cut_i];
                        double r = std::ceil(sum_demands / C);
                        double current_cut_value = 0.0;
                        GRBLinExpr c = 0.0;

                        for (auto i = cuts[cut_i].begin(); i != cuts[cut_i].end(); ++i)
                            for (auto j = std::next(i); j != cuts[cut_i].end(); ++j)
                            {
                                c += x[std::max(*i,*j)][std::min(*i,*j)];
                                current_cut_value += getNodeRel(x[std::max(*i,*j)][std::min(*i,*j)]);
                            }

                        if (current_cut_value > cuts[cut_i].size() - r)
                            addLazy(c, GRB_LESS_EQUAL, cuts[cut_i].size() - r);
                    }
                }

                cut_generator.clear_edges();
            }
        }
};

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: vrp <filename>" << std::endl;
        return 1;
    }
    int V, N;
    double C;

    std::vector<Point> clients = getPointsFromFile(argv[1], &V, &C);
    N = clients.size();
    double **heur_vals = getHeuristicSol(argv[1], N, V);

    try
    {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);
        model.set(GRB_DoubleParam_TimeLimit, 1 * 60 * 60);

        // x_e variables
        GRBVar **x = new GRBVar*[N];
        x[0] = nullptr;
        for (int i = 1; i < N; i++)
        {
            GRBVar *row = new GRBVar[i];
            for (int j = 0; j < i; j++)
            {
                std::string varname = "x[" + std::to_string(i) + "][" + std::to_string(j) + "]";
                double upper_limit = j == 0 ? 2.0 : 1.0;
                row[j] = model.addVar(0.0, upper_limit, dist(clients[i], clients[j]), GRB_INTEGER, varname);
                row[j].set(GRB_DoubleAttr_Start, heur_vals[i][j]);
            }
            x[i] = row;
        }

        // Free heuristic values space
        for (int i = 1; i < N; i++)
            delete[] heur_vals[i];
        delete[] heur_vals;

        // Degree constraints (clients)
        for (int i = 1; i < N; i++)
        {
            GRBLinExpr c = 0.0;
            for (int j = 0; j < N; j++)
                if (i != j)
                    c += x[std::max(i,j)][std::min(i,j)];
            model.addConstr(c, GRB_EQUAL, 2);
        }

        // Degree constraints (warehouse)
        GRBLinExpr c = 0.0;
        for (int i = 1; i < N; i++)
            c += x[i][0];
        model.addConstr(c, GRB_LESS_EQUAL, 2 * V);

        double* demands = new double[N];
        for (int i = 0; i < N; i++)
            demands[i] = clients[i].d;

        subtourelim cb(x, demands, N, V, C);
        model.setCallback(&cb);
        model.optimize();

        printf("\nGap: %lf\nRuntime: %lf\n", model.get(GRB_DoubleAttr_MIPGap), model.get(GRB_DoubleAttr_Runtime));
        writeSolution(x, N, V, model.get(GRB_DoubleAttr_ObjVal), 1, argv[1]);

        // Deallocating
        for (int i = 1; i < N; i++)
            delete[] x[i];
        delete[] x;
        delete[] demands;
    } 
    catch (GRBException e) 
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } 
    catch (...) 
    {
        std::cout << "Unknown error" << std::endl;
    }

    return 0;
}

double dist(const Point& a, const Point& b)
{
    return std::sqrt( (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) );
}

std::vector<Point> getPointsFromFile(std::string filename, int *V, double *C)
{
    int N;
    
    std::ifstream f(filename);
    std::string linebuffer;
    if (!f.is_open()) 
    {
        std::cout << "Error: couldn't open file " << filename << std::endl;
        exit(1);
    }

    getline(f, linebuffer);
    sscanf(linebuffer.c_str(), "%i %i %lf\n", &N, V, C);

    std::vector<Point> points;
    points.reserve(N);

    while (getline(f, linebuffer)) 
    {
        if (!linebuffer.empty()) 
        {
            Point p;
            sscanf(linebuffer.c_str(), "%lf %lf %lf", &p.d, &p.x, &p.y);
            points.push_back(p);
        }
    }

    f.close();
    return points;
}

double** getHeuristicSol(std::string filename, int N, int V)
{
    std::ifstream f(filename + ".heu");
    std::string linebuffer;
    if (!f.is_open())
    {
        std::cout << "Error: couldn't open file " << filename << std::endl;
        exit(1);
    }

    double** vals = new double*[N];
    vals[0] = nullptr;
    for (int i = 1; i < N; i++)
    {
        vals[i] = new double[i];
        for (int j = 0; j < i; j++)
            vals[i][j] = 0.0;
    }

    getline(f, linebuffer);
    for (int k = 0; k < V; k++)
    {
        int last_vertex, current_vertex;
        f >> last_vertex;
        f >> current_vertex;
        while (current_vertex != 0)
        {
            vals[std::max(last_vertex, current_vertex)][std::min(last_vertex, current_vertex)] += 1.0;
            last_vertex = current_vertex;
            f >> current_vertex;
        }

        if (last_vertex != current_vertex)
            vals[std::max(last_vertex, current_vertex)][std::min(last_vertex, current_vertex)] += 1.0;
    }

    f.close();
    return vals;
}

void writeSolution(GRBVar** x, int N, int V, const double obj, const int opt, const std::string& filename)
{
    std::ofstream f(filename + ".sol");
    if (!f.is_open()) 
    {
        std::cout << "Error: couldn't open file " << filename + ".sol" << std::endl;
        exit(1);
    }

    f << std::fixed << std::setprecision(2) << obj << ' ' << opt << '\n';
    std::vector<bool> seen;
    for (int i = 0; i < N; i++)
        seen.push_back(false);

    for (int k = 0; k < V; k++)
    {
        std::list<int> tour;
        tour.push_back(0);
        int next;
        int current = 0;
        do {
            next = 0;
            for (int j = 1; j < N; j++)
            {
                if (current != j && !seen[j] && x[std::max(current, j)][std::min(current, j)].get(GRB_DoubleAttr_X) > 0.5)
                {
                    next = j;
                    tour.push_back(j);
                    current = j;
                    seen[j] = true;
                    break;
                }
            }
        } while (next != 0);
        tour.push_back(0);

        for (auto it = tour.begin(); it != tour.end(); ++it)
        {
            f << *it;
            if (std::next(it) == tour.end())
                f << '\n';
            else
                f << ' ';
        }
    }

    f.close();
}