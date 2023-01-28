# Super basic makefile for examples.

CXX=g++ -std=c++17
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion
OPT=-O3

all: example Bessel

clean:
	rm -f -v example

EXDIR=./examples

example: $(EXDIR)/example.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I. -lgsl

Bessel: $(EXDIR)/Bessel.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I. -lgsl