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


all: reduce_rand reduce_band

reduce_rand: $(OBJS)
	$(MPICXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_rand.cpp -o reduce_rand

reduce_band: $(OBJS)
	$(MPICXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_band.cpp -o reduce_band


clean:
	rm -rf reduce_rand reduce_band
