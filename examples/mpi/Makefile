-include ../make.inc
-include ../common.inc
MPICXX = mpicxx
OMPI_CXX = $(CXX)
CPPFLAGS = -I../../src


all: clean spmat

spmat: $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spmat.cpp -o spmat

clean:
	rm -rf spmat
