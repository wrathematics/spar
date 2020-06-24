-include make.inc
CPPFLAGS = -I../src
#DEBUGFLAGS = -g


all: reduce_rand reduce_bandish

reduce_rand: $(OBJS)
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_rand.cpp -o reduce_rand

reduce_bandish: $(OBJS)
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_bandish.cpp -o reduce_bandish


clean:
	rm -rf reduce_rand reduce_bandish
