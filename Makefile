GRBPATH = /opt/gurobi902/linux64
LIBS = -lgurobi_c++ -lgurobi90
FLAGS = -g -Wall -pedantic -Wno-unused-result

unionfind.o: unionfind.cpp unionfind.hpp
	g++ $(FLAGS) -c unionfind.cpp -o unionfind.o

karger.o: karger.cpp karger.hpp
	g++ $(FLAGS) -c karger.cpp -o karger.o

kruskal.o: kruskal.cpp kruskal.hpp
	g++ $(FLAGS) -c kruskal.cpp -o kruskal.o

exec: vrp.cpp karger.o unionfind.o kruskal.o
	g++ $(FLAGS) vrp.cpp karger.o unionfind.o kruskal.o -I$(GRBPATH)/include -L$(GRBPATH)/lib $(LIBS) -o CVRPSolver