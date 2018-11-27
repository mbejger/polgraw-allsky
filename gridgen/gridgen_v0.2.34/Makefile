CXX := g++
INCDIR := -I./include 
CXXFLAGS := -g -std=c++11 -Wall -Wextra -O2 -s

SRCFILES := $(wildcard main.cpp src/*.cpp) 
OBJFILES := $(patsubst %.cpp,%.o,$(SRCFILES))

gridgen: $(OBJFILES)
	$(CXX) $^ -o $@

%.o: %.cpp
	$(CXX) $(INCDIR) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o src/*.o gridgen
