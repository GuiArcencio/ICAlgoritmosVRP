GRBPATH = /opt/gurobi902/linux64
LIBS = -lgurobi_c++ -lgurobi90
FLAGS = -O3 -g -Wall -pedantic -Wno-unused-result

unionfind.o: unionfind.cpp unionfind.hpp
	g++ $(FLAGS) -c unionfind.cpp -o unionfind.o

karger.o: karger.cpp karger.hpp
	g++ $(FLAGS) -c karger.cpp -o karger.o

exec: vrp.cpp karger.o unionfind.o
	g++ $(FLAGS) vrp.cpp karger.o unionfind.o -I$(GRBPATH)/include -L$(GRBPATH)/lib $(LIBS) -o CVRPSolver