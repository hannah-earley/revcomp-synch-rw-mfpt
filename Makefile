ifeq (, $(shell which clang++ 2>/dev/null))
    CXX=g++
    CXXFLAGS=-fopenmp -std=c++11 -O3 -Wall -Iinc/
    LINKFLAGS=
else
    CXX=clang++
    CXXFLAGS=-Xclang -fopenmp -std=c++11 -O3 -Wall -Iinc/
    LINKFLAGS=-lomp
endif
SRCDIR=src
OBJDIR=obj

_DEPS=walk.h testbed.cpp
DEPS=$(patsubst %,$(SRCDIR)/%,$(_DEPS))
_OBJ=walk.o misc.o
OBJ=$(patsubst %,$(OBJDIR)/%,$(_OBJ))

all: walk

walk: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: walk
	./walk

clean:
	-rm -f $(OBJDIR)/*.o
	-rmdir $(OBJDIR)


.PHONY: all run clean
