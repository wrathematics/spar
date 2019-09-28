-include make.inc

CPPFLAGS = -I../src
CXXFLAGS = -std=c++11 -O2 -march=native
OMPFLAGS = -fopenmp 
WARNFLAGS = -Wall -pedantic -Wextra


all: miniclean spvec run

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) -c $< -o $@

spvec: $(OBJS)
	$(CXX) $(OBJS) -o spvec $(OMPFLAGS)

run: spvec
	( ./spvec )

miniclean:
	rm -rf spvec basics.o

clean:
	rm -rf *.o spvec
