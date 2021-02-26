#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <list>
#include <limits>

#include "gurobi_c++.h"
#include "karger.hpp"
#include "kruskal.hpp"

#include "cxxopts/cxxopts.hpp"

struct Point 
{
    double x;
    double y;
    double d;
};

double dist(const Point& a, const Point& b);

std::vector<Point> getPointsFromFile(std::string filename, int *V, double *C);
double** getHeuristicSol(std::string filename, int N, int V, double* upper_bound);
void writeSolution(GRBVar** x, int N, int V, const double obj, const int opt, const std::string& filename);

class subtourelim: public GRBCallback
{
    public:
        GRBVar** x;
        int N, V;
        double C, coefficient;
        double* demands;
        bool use_log;
        karger::EdgeVector cut_generator;
        subtourelim(GRBVar** x, double* demands, int N, int V, double C, double coefficient, bool use_log):
            x(x), N(N), V(V), C(C), demands(demands), cut_generator(N, 1), coefficient(coefficient), use_log(use_log) {};
        

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
                std::vector<kruskal::edge> kruskal_edges;

                // Setup
                for (int i = 1; i < N; i++)
                    for (int j = 0; j < i; j++)
                    {
                        double rel_val = getNodeRel(x[i][j]);

                        // Karger
                        if (j > 0)
                            cut_generator.add_edge(i, j, rel_val);

                        // Kruskal
                        kruskal::edge e;
                        e.u = i;
                        e.v = j;
                        e.length = x[i][j].get(GRB_DoubleAttr_Obj);
                        e.starting = rel_val > 0.9;

                        kruskal_edges.push_back(e);
                    }

                // Run Karger 
                for (int count = 0; count < coefficient * (use_log ? std::log(N) : N); count++)
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

                // --------------- GERAÇÃO DE COVER CUTS --------------------------------
                
                double msg_cost, max_length;
                auto msg = kruskal::generateMSG(N, V, kruskal_edges, &msg_cost, &max_length);

                if (msg_cost > getDoubleInfo(GRB_CB_MIPNODE_OBJBND))
                {
                    bool** in_msg = new bool*[N];
                    in_msg[0] = nullptr;
                    for (int i = 1; i < N; i++)
                    {
                        in_msg[i] = new bool[i];
                        for (int j = 0; j < i; j++)
                            in_msg[i][j] = false;
                    }

                    GRBLinExpr c = 0.0;
                    for (auto e : msg)
                    {
                        c += x[e.u][e.v];
                        in_msg[e.u][e.v] = true;
                    }

                    // Extended cover
                    for (auto e : kruskal_edges)
                        if (!in_msg[e.u][e.v] && e.length > max_length)
                            c += x[e.u][e.v];

                    addLazy(c, GRB_LESS_EQUAL, msg.size() - 1);

                    for (int i = 1; i < N; i++)
                        delete[] in_msg[i];
                    delete[] in_msg;
                }
            }
        }
};

int main(int argc, char *argv[])
{
    int V, N;
    double C, time_limit, upper_bound = std::numeric_limits<double>::infinity();
    double coefficient;
    bool use_log, use_heur;
    std::vector<Point> clients;
    double **heur_vals = nullptr;

    cxxopts::Options options("CVRPSolver", "Outputs an exact solution for a CVRP instance");
    options.add_options()
        ("f,file", "Input file name", cxxopts::value<std::string>())
        ("H,use-heuristic", "Use a heuristic solution from <input-file-name>.heu", cxxopts::value<bool>()->default_value("false"))
        ("C,karger-coefficient", "Constant coefficient for how many times Karger's Algorithm will be executed", cxxopts::value<double>()->default_value("10.0"))
        ("l,use-log-n", "Run Karger's Algorithm O(log n) times instead of O(n)", cxxopts::value<bool>()->default_value("false"))
        ("t,time-limit", "Time limit for the solver in seconds", cxxopts::value<double>()->default_value("3600.0"))
        ("h,help", "Prints this page")
    ;

    auto command_line = options.parse(argc, argv);

    try 
    {

    if (command_line.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (!command_line.count("file"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    coefficient = command_line["karger-coefficient"].as<double>();
    use_log = command_line["use-log-n"].as<bool>();
    use_heur = command_line["use-heuristic"].as<bool>();
    time_limit = command_line["time-limit"].as<double>();

    clients = getPointsFromFile(command_line["file"].as<std::string>(), &V, &C);
    N = clients.size();

    if (use_heur)
        heur_vals = getHeuristicSol(command_line["file"].as<std::string>(), N, V, &upper_bound);

    }
    catch(cxxopts::OptionException* e) {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    try
    {
        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);
        model.set(GRB_DoubleParam_TimeLimit, time_limit);

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
                if (use_heur)
                    row[j].set(GRB_DoubleAttr_Start, heur_vals[i][j]);
            }
            x[i] = row;
        }

        // Free heuristic values space
        if (use_heur)
        {
            for (int i = 1; i < N; i++)
                delete[] heur_vals[i];
            delete[] heur_vals;
        }

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

        subtourelim cb(x, demands, N, V, C, coefficient, use_log);
        model.setCallback(&cb);
        model.optimize();

        printf("\nOBJ: %lf\nGap: %lf\nRuntime: %lf\n", model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_MIPGap), model.get(GRB_DoubleAttr_Runtime));
        writeSolution(x, N, V, model.get(GRB_DoubleAttr_ObjVal), 1, command_line["file"].as<std::string>());

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

double** getHeuristicSol(std::string filename, int N, int V, double* upper_bound)
{
    std::ifstream f(filename + ".heu");
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

    f >> (*upper_bound);
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