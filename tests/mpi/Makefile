-include ../make.inc
-include ../common.inc
MPICXX = mpicxx
OMPI_CXX = $(CXX)
CPPFLAGS = -I../../src -I../


all: mpi

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

%.o: %.cpp 
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) -c $< -o $@
 
mpi: $(OBJS)
	$(MPICXX) $(OBJS) -o mpi $(OMPFLAGS)

run: mpi
	( mpirun -np 2 ./mpi )

miniclean:
	( ls *.o | grep -v test-main | xargs rm -f )

clean:
	rm -rf *.o mpi
