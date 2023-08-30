CXX       ?= g++
CXXFLAGS  ?= -std=c++17
LINK.o    := $(LINK.cc) # Use C++ linker

CPPFLAGS  += -fopenmp -O3 -Wall -Wextra -I. -I../Mesh -I../VirtualSpace -I${mkBoostInc} -I${mkEigenInc}
LDFLAGS = -L. -L${mkBoostLib} -L${mkEigenLib} -Wl,-rpath=${PWD}
LIBS += -lboost_system -lboost_filesystem

# Source files and object files for main_integration
SOLVER_SRCS  = main_solver.cpp monomial.cpp
SOLVER_OBJS  = $(SOLVER_SRCS:.cpp=.o)

# Source files and object files for main_exportresults
EXPORTRESULTS_SRCS  = main_exportresults.cpp
EXPORTRESULTS_OBJS  = $(EXPORTRESULTS_SRCS:.cpp=.o)

# Source files and object files for main_parameters
PARAMETERS_SRCS  = main_parameters.cpp
PARAMETERS_OBJS  = $(PARAMETERS_SRCS:.cpp=.o)

# Source files and object files for main
MAIN_SRCS  = main.cpp monomial.cpp
MAIN_OBJS  = $(MAIN_SRCS:.cpp=.o)

HEADERS = $(wildcard *.hpp *.h *.hh)

# Executable names
SOLVER_EXEC     	:= main_solver
EXPORTRESULTS_EXEC	:= main_exportresults
PARAMETERS_EXEC 	:= main_parameters
MAIN_EXEC 			:= main

.PHONY: all clean main_solver main_exportresults.cpp main_parameters.cpp main.cpp

.DEFAULT_GOAL := all

# Targets and rules
all: $(SOLVER_EXEC) $(EXPORTRESULTS_EXEC) $(PARAMETERS_EXEC) $(MAIN_EXEC)

# Targets and rules for main_solver
$(SOLVER_EXEC): $(SOLVER_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main_exportresults
$(EXPORTRESULTS_EXEC): $(EXPORTRESULTS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main_parameters
$(PARAMETERS_EXEC): $(PARAMETERS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main
$(MAIN_EXEC): $(MAIN_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Compile .cpp files into .o files
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

# Clean up object files and executables
clean:
	rm -f $(SOLVER_OBJS) $(SOLVER_EXEC) $(EXPORTRESULTS_OBJS) $(EXPORTRESULTS_EXEC) \
		  $(PARAMETERS_OBJS) $(PARAMETERS_EXEC) $(MAIN_OBJS) $(MAIN_EXEC)
