# Super basic makefile for examples.
CXX=g++ 
STD=-std=c++17
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion
OPT=-O3

all: Bessel Complex Schrodinger Dirac Inhomogenous

clean:
	rm -f -v Bessel Complex Schrodinger Dirac Inhomogenous test

EXDIR=./examples

# This one uses GSL library, to check against exact Bessel solution
Bessel: $(EXDIR)/Bessel.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I. -lgsl

Complex: $(EXDIR)/Complex.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.

Schrodinger: $(EXDIR)/Schrodinger.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.

Dirac: $(EXDIR)/Dirac.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.

Inhomogenous: $(EXDIR)/Inhomogenous.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ $(EXDIR)/$@.cpp $(WARN) $(OPT) -I.

# Build the unit test executable (not built by default)
test: ./tests/test.cpp AdamsMoulton.hpp
	$(CXX) $(STD) -o $@ ./tests/$@.cpp $(WARN) $(OPT) -I.

# performs code coverage (test) report, for code-cov
.PHONY: coverage
coverage:
	g++-11 --coverage -g -o coverage ./tests/test.cpp -I.
	./coverage
	lcov --capture --directory . --output-file coverage.info
	lcov --remove coverage.info '*/catch2/*' '*/tests/*' '*/examples/*' '/usr/*' --output-file coverage.info
	lcov --list coverage.info |tee -a cov-info.txt