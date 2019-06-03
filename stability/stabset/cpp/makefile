CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs   ibex)

CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-unused-variable -fopenmp -frounding-math -std=c++11 -I/home/alexandre/Eigen -I/home/alexandre/libigl/include/igl
EXEC= test_optim_kd

all: $(EXEC)

test_optim_kd: kharitonov.o dabbene.o optim_KD.o test_optim_kd.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	g++ -o $@ -c $< $(CXXFLAGS) $(LIBS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
