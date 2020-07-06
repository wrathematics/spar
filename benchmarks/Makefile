MPICXX = mpicxx
OMPI_CXX = g++

CXXFLAGS = -std=c++17 -O2
OMPFLAGS = -fopenmp 
WARNFLAGS = -Wall -pedantic -Wextra

CPPFLAGS = -I../src/spar
#DEBUGFLAGS = -g

ifneq ($(ARCH),$(filter $(ARCH),ppc64 ppc64le))
CXXFLAGS += -mcpu=native
else
CXXFLAGS += -march=native
endif


all: reduce_rand reduce_bandish

reduce_rand: $(OBJS)
	$(MPICXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_rand.cpp -o reduce_rand

reduce_bandish: $(OBJS)
	$(MPICXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_bandish.cpp -o reduce_bandish


clean:
	rm -rf reduce_rand reduce_bandish
