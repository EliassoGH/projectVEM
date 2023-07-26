CXX      ?= g++
CXXFLAGS ?= -std=c++20
LINK.o := $(LINK.cc) # Use C++ linker

CPPFLAGS += -O3 -Wall -I.

# Source files and object files
POINT_SRCS := main_point.cpp
EDGE_SRCS := main_edge.cpp
POLYGON_SRCS := main_polygon.cpp
POLYHEDRON_SRCS := main_polyhedron.cpp
MESH_SRCS := main_mesh.cpp

POINT_OBJS := $(POINT_SRCS:.cpp=.o)
EDGE_OBJS := $(EDGE_SRCS:.cpp=.o)
POLYGON_OBJS := $(POLYGON_SRCS:.cpp=.o)
POLYHEDRON_OBJS := $(POLYHEDRON_SRCS:.cpp=.o)
MESH_OBJS := $(MESH_SRCS:.cpp=.o)

HEADERS = $(wildcard *.hpp *.h *.hh)

# Executable names
POINT_EXEC := main_point
EDGE_EXEC := main_edge
POLYGON_EXEC := main_polygon
POLYHEDRON_EXEC := main_polyhedron
MESH_EXEC := main_mesh

.PHONY = all clean main_point main_edge main_polygon main_polyhedron main_mesh

.DEFAULT_GOAL = all

# Targets and rules
all: $(POINT_EXEC) $(EDGE_EXEC) $(POLYGON_EXEC) $(POLYHEDRON_EXEC) $(MESH_EXEC)

$(POINT_EXEC): $(POINT_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EDGE_EXEC): $(EDGE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(POLYGON_EXEC): $(POLYGON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(POLYHEDRON_EXEC): $(POLYHEDRON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(MESH_EXEC): $(MESH_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile .cpp files into .o files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up object files and executables
clean:
	rm -f $(POINT_OBJS) $(EDGE_OBJS) $(POLYGON_OBJS) $(POLYHEDRON_OBJS) $(MESH_OBJS) \
		  $(POINT_EXEC) $(EDGE_EXEC) $(POLYGON_EXEC) $(POLYHEDRON_EXEC) $(MESH_EXEC)
