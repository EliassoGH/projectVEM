CXX      ?= g++
CXXFLAGS ?= -std=c++20
LINK.o := $(LINK.cc) # Use C++ linker

CPPFLAGS += -O3 -Wall -I.

# Source files and object files
POINT_SRCS := main_point.cpp
EDGE_SRCS := main_edge.cpp

POINT_OBJS := $(POINT_SRCS:.cpp=.o)
EDGE_OBJS := $(EDGE_SRCS:.cpp=.o)

HEADERS = $(wildcard *.hpp *.h *.hh)

# Executable names
POINT_EXEC := main_point
EDGE_EXEC := main_edge

.PHONY = all clean

.DEFAULT_GOAL = all

# Targets and rules
all: $(POINT_EXEC) $(EDGE_EXEC)

$(POINT_EXEC): $(POINT_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EDGE_EXEC): $(EDGE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile .cpp files into .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up object files and executables
clean:
	rm -f $(POINT_OBJS) $(EDGE_OBJS) $(POINT_EXEC) $(EDGE_EXEC)
