-include ../make.inc
-include ../common.inc

MPICXX = mpicxx
OMPI_CXX = $(CXX)
EIGEN_CPPFLAGS = -I/usr/include/eigen3/


all: clean eigen spmat

eigen: $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(EIGEN_CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) eigen.cpp -o eigen

spmat: $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spmat.cpp -o spmat

clean:
	rm -rf eigen spmat
