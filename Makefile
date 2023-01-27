# Super basic makefile for examples.

CXX=g++ -std=c++17
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion
OPT=-O3

all: example

clean:
	rm -f -v example

example: example.cpp AdamsMoulton.hpp
	$(CXX) -o $@ $@.cpp $(WARN) $(OPT) -lgsl

