-include ../make.inc
-include ../common.inc


all: clean gen spmat spvec

gen: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) gen.cpp -o gen

spmat: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spmat.cpp -o spmat

spvec: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) spvec.cpp -o spvec

clean:
	rm -rf gen spmat spvec
