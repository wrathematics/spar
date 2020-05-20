-include make.inc
CPPFLAGS = -I../src


all: clean reduce_rand

reduce_rand: $(OBJS)
	$(CXX) -g $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) reduce_rand.cpp -o reduce_rand

clean:
	rm -rf reduce_rand
