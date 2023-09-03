CXX      ?= g++
CXXFLAGS ?= -std=c++17 -g
#LINK.o := $(LINK.cc) # Use C++ linker
CPPFLAGS += -fopenmp -O3 -Wall -Iinclude -I${mkBoostInc} -I${mkEigenInc}

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files and object files
POINT_SRCS := $(SRC_DIR)/main_point.cpp
EDGE_SRCS := $(SRC_DIR)/main_edge.cpp
POLYGON_SRCS := $(SRC_DIR)/main_polygon.cpp
POLYHEDRON_SRCS := $(SRC_DIR)/main_polyhedron.cpp
MESH_SRCS := $(SRC_DIR)/main_mesh.cpp

MONOMIAL_SRCS := $(SRC_DIR)/main_monomial.cpp $(SRC_DIR)/monomial.cpp
INTEGRATION_SRCS := $(SRC_DIR)/main_integration.cpp $(SRC_DIR)/monomial.cpp $(SRC_DIR)/integration.cpp
VIRTUALDOFS_SRCS := $(SRC_DIR)/main_virtualdofs.cpp $(SRC_DIR)/monomial.cpp \
					$(SRC_DIR)/integration.cpp $(SRC_DIR)/virtualDofs.cpp
VIRTUALPROJECTIONS_SRCS := $(SRC_DIR)/main_virtualprojections.cpp $(SRC_DIR)/monomial.cpp \
						   $(SRC_DIR)/integration.cpp $(SRC_DIR)/virtualDofs.cpp \
						   $(SRC_DIR)/virtualProjections.cpp

PARAMETERS_SRCS := $(SRC_DIR)/main_parameters.cpp $(SRC_DIR)/parameters.cpp
MAIN_SRCS := $(SRC_DIR)/main.cpp $(SRC_DIR)/monomial.cpp $(SRC_DIR)/integration.cpp \
			 $(SRC_DIR)/virtualDofs.cpp $(SRC_DIR)/virtualProjections.cpp \
			 $(SRC_DIR)/solver.cpp $(SRC_DIR)/problem.cpp $(SRC_DIR)/export_results.cpp \
			 $(SRC_DIR)/parameters.cpp


POINT_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(POINT_SRCS))
EDGE_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(EDGE_SRCS))
POLYGON_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(POLYGON_SRCS))
POLYHEDRON_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(POLYHEDRON_SRCS))
MESH_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(MESH_SRCS))

MONOMIAL_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(MONOMIAL_SRCS))
INTEGRATION_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(INTEGRATION_SRCS))
VIRTUALDOFS_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(VIRTUALDOFS_SRCS))
VIRTUALPROJECTIONS_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(VIRTUALPROJECTIONS_SRCS))

PARAMETERS_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(PARAMETERS_SRCS))
MAIN_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(MAIN_SRCS))

HEADERS = $(wildcard *.hpp *.h *.hh)

# Executable names
POINT_EXEC := main_point
EDGE_EXEC := main_edge
POLYGON_EXEC := main_polygon
POLYHEDRON_EXEC := main_polyhedron
MESH_EXEC := main_mesh

MONOMIAL_EXEC := main_monomial
INTEGRATION_EXEC := main_integration
VIRTUALDOFS_EXEC := main_virtualdofs
VIRTUALPROJECTIONS_EXEC := main_virtualprojections

PARAMETERS_EXEC := main_parameters
MAIN_EXEC := main


.PHONY = all clean main_point main_edge main_polygon main_polyhedron main_mesh \
		main_monomial main_integration main_virtualdofs main_virtualprojections \
		main_parameters main

.DEFAULT_GOAL = main

# Targets and rules
all: $(POINT_EXEC) $(EDGE_EXEC) $(POLYGON_EXEC) $(POLYHEDRON_EXEC) $(MESH_EXEC) \
	 $(MONOMIAL_EXEC) $(INTEGRATION_EXEC) $(VIRTUALDOFS_EXEC) $(VIRTUALPROJECTIONS_EXEC) \
	 $(PARAMETERS_EXEC) $(MAIN_EXEC)

$(POINT_EXEC): $(POINT_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(EDGE_EXEC): $(EDGE_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(POLYGON_EXEC): $(POLYGON_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(POLYHEDRON_EXEC): $(POLYHEDRON_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(MESH_EXEC): $(MESH_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^


$(MONOMIAL_EXEC): $(MONOMIAL_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(INTEGRATION_EXEC): $(INTEGRATION_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(VIRTUALDOFS_EXEC): $(VIRTUALDOFS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(VIRTUALPROJECTIONS_EXEC): $(VIRTUALPROJECTIONS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^


$(PARAMETERS_EXEC): $(PARAMETERS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(MAIN_EXEC): $(MAIN_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^


# Compile .cpp files into .o files
#%.o: %.cpp
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

# Clean up object files and executables
clean:
	rm -f $(POINT_OBJS) $(EDGE_OBJS) $(POLYGON_OBJS) $(POLYHEDRON_OBJS) $(MESH_OBJS) \
		  $(POINT_EXEC) $(EDGE_EXEC) $(POLYGON_EXEC) $(POLYHEDRON_EXEC) $(MESH_EXEC) \
		  $(MONOMIAL_EXEC) $(INTEGRATION_EXEC) $(VIRTUALDOFS_EXEC) $(VIRTUALPROJECTIONS_EXEC) \
		  $(MONOMIAL_OBJS) $(INTEGRATION_OBJS) $(VIRTUALDOFS_OBJS) $(VIRTUALPROJECTIONS_OBJS) \
		  $(PARAMETERS_OBJS) $(MAIN_OBJS)\
		  $(PARAMETERS_EXEC) $(MAIN_EXEC)
