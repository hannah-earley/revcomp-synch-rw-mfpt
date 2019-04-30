ifeq (, $(shell which clang++))
    CXX=g++
    CXXFLAGS=-fopenmp -std=c++11 -O3 -Wall -Iinc/
    LINKFLAGS=
else
    CXX=clang++
    CXXFLAGS=-Xclang -fopenmp -std=c++11 -O3 -Wall -Iinc/
    LINKFLAGS=-lomp
endif
OBDIR=obj

DEPS=walk.h testbed.cpp
_OBJ=walk.o misc.o
OBJ=$(patsubst %,$(OBDIR)/%,$(_OBJ))

all: walk

walk: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) -o $@ $^

$(OBDIR)/%.o: %.cpp $(DEPS)
	mkdir -p $(OBDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: walk
	./walk

clean:
	-rm -f $(OBDIR)/*.o
	-rmdir $(OBDIR)


.PHONY: all run clean
