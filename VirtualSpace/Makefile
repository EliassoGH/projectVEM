CXX       ?= g++
CXXFLAGS  ?= -std=c++17
LINK.o    := $(LINK.cc) # Use C++ linker

CPPFLAGS  += -fopenmp -O3 -Wall -Wextra -I. -I../Mesh -I../Solver -I${mkBoostInc} -I${mkEigenInc}
LDFLAGS = -L. -L${mkBoostLib} -L${mkEigenLib} -Wl,-rpath=${PWD}
LIBS += -lboost_system -lboost_filesystem

# Source files and object files for main_monomial
MONOMIAL_SRCS  = main_monomial.cpp
MONOMIAL_OBJS  = $(MONOMIAL_SRCS:.cpp=.o)

# Source files and object files for main_virtualdofs
VIRTUALDOFS_SRCS  = main_virtualdofs.cpp
VIRTUALDOFS_OBJS  = $(VIRTUALDOFS_SRCS:.cpp=.o)

# Source files and object files for main_integration
INTEGRATION_SRCS  = main_integration.cpp
INTEGRATION_OBJS  = $(INTEGRATION_SRCS:.cpp=.o)

# Source files and object files for main_integration
VIRTUALPROJECTIONS_SRCS  = main_virtualprojections.cpp
VIRTUALPROJECTIONS_OBJS  = $(VIRTUALPROJECTIONS_SRCS:.cpp=.o)

HEADERS = $(wildcard *.hpp *.h *.hh)

# Executable names
MONOMIAL_EXEC     := main_monomial
VIRTUALDOFS_EXEC  := main_virtualdofs
INTEGRATION_EXEC  := main_integration
VIRTUALPROJECTIONS_EXEC  := main_virtualprojections

.PHONY: all clean main_monomial main_virtualdofs main_integration main_virtualprojections

.DEFAULT_GOAL := all

# Targets and rules
all: $(MONOMIAL_EXEC) $(VIRTUALDOFS_EXEC) $(INTEGRATION_EXEC) $(VIRTUALPROJECTIONS_EXEC)

# Targets and rules for main_monomial
$(MONOMIAL_EXEC): $(MONOMIAL_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main_virtualdofs
$(VIRTUALDOFS_EXEC): $(VIRTUALDOFS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main_integration
$(INTEGRATION_EXEC): $(INTEGRATION_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Targets and rules for main_virtualprojections
$(VIRTUALPROJECTIONS_EXEC): $(VIRTUALPROJECTIONS_OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Compile .cpp files into .o files
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

# Clean up object files and executables
clean:
	rm -f $(MONOMIAL_OBJS) $(MONOMIAL_EXEC) $(VIRTUALDOFS_OBJS) $(VIRTUALDOFS_EXEC) \
		  $(INTEGRATION_OBJS) $(INTEGRATION_EXEC) $(VIRTUALPROJECTIONS_OBJS) $(VIRTUALPROJECTIONS_EXEC)
