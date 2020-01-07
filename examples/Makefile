-include make.inc

CPPFLAGS = -I../src
CXXFLAGS = -std=c++17 -O2 -march=native
OMPFLAGS = -fopenmp 
WARNFLAGS = -Wall -pedantic -Wextra


all: clean gen spmat spvec

gen: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) gen.cpp -o gen

spmat: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spmat.cpp -o spmat

spvec: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spvec.cpp -o spvec

clean:
	rm -rf gen spmat spvec
