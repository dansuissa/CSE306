CXX = g++
CXXFLAGS = -O3 -std=c++11 -fopenmp -Wall -Wextra -I.
LDFLAGS = -fopenmp
all: voronoi optimal_transport

voronoi: cvt.o Voronoi.o
	$(CXX) $(CXXFLAGS) -o voronoi cvt.o Voronoi.o $(LDFLAGS)

optimal_transport: ot_gaussian.o PowerDiagram.o lbfgs.o
	$(CXX) $(CXXFLAGS) -o optimal_transport ot_gaussian.o PowerDiagram.o lbfgs.o $(LDFLAGS)

cvt.o: cvt.cpp
	$(CXX) $(CXXFLAGS) -c cvt.cpp

ot_gaussian.o: ot_gaussian.cpp
	$(CXX) $(CXXFLAGS) -c ot_gaussian.cpp

Voronoi.o: Voronoi.cpp
	$(CXX) $(CXXFLAGS) -c Voronoi.cpp

PowerDiagram.o: PowerDiagram.cpp
	$(CXX) $(CXXFLAGS) -c PowerDiagram.cpp

lbfgs.o: lbfgs.c
	$(CXX) $(CXXFLAGS) -c lbfgs.c

clean:
	rm -f *.o voronoi optimal_transport 