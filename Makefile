CXX=clang++
CXXFLAGS=-std=c++11 -O3 -Iinc/
OBDIR=obj

DEPS=
_OBJ=walk.o
OBJ=$(patsubst %,$(OBDIR)/%,$(_OBJ))

all: walk

walk: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBDIR)/%.o: %.cpp $(DEPS)
	mkdir -p $(OBDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: walk
	./walk

clean:
	rm -f $(OBDIR)/*.o
	rmdir $(OBDIR)


.PHONY: all run clean