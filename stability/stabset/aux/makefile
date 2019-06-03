CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs   ibex)

CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-unused-variable -fopenmp -frounding-math -std=c++11 -I/home/alexandre/Eigen -I/home/alexandre/libigl/include/igl -std=c++11

EXEC= test_optim_kd

all: $(EXEC)

test_topology: kharitonov.o test_topology.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

test_stability: kharitonov.o dabbene.o test_stability.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

test_kharitonov: kharitonov.o test_kharitonov.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)
	
test_dabbene: dabbene.o test_dabbene.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	g++ -o $@ -c $< $(CXXFLAGS) $(LIBS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
