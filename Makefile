# Super basic makefile for examples.

CXX=g++ -std=c++17
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion
OPT=-O3

all: Bessel Complex Schrodinger

clean:
	rm -f -v example Bessel Complex Schrodinger

EXDIR=./examples

# This one uses GSL library, to check against exact Bessel solution
Bessel: $(EXDIR)/Bessel.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I. -lgsl

Complex: $(EXDIR)/Complex.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.

Schrodinger: $(EXDIR)/Schrodinger.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.