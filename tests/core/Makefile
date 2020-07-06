-include ../make.inc
-include ../common.inc


all: core

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(WARNFLAGS) $(OMPFLAGS) -c $< -o $@
 
core: $(OBJS)
	$(CXX) $(OBJS) -o core $(OMPFLAGS)

run: core
	( ./core )

miniclean:
	( ls *.o | grep -v test-main | xargs rm -f )

clean:
	rm -rf *.o core
