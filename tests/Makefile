-include make.inc

CPPFLAGS = -I../src
CXXFLAGS = -std=c++17 -O2 -march=native
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
	( ls *.o | grep -v test-main | xargs rm -f )

clean:
	rm -rf *.o spvec
