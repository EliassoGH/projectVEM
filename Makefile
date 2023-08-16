CXX      ?= g++
CXXFLAGS ?= -std=c++20 -g
LINK.o := $(LINK.cc) # Use C++ linker

CPPFLAGS += -O3 -Wall -I. -I./GTMesh/src/

SRCS = $(wildcard *.cpp *.cc)
OBJS = $(SRCS:.cpp=.o)
HEADERS = $(wildcard *.hpp *.h *.hh)

exe_sources = $(filter main%.cpp,$(SRCS))
EXEC = $(exe_sources:.cpp=)

.PHONY = all parallel clean distclean

.DEFAULT_GOAL = all

all: $(EXEC)

$(EXEC): $(OBJS)

$(OBJS): $(SRCS) $(HEADERS)

clean:
	$(RM) -f $(EXEC) $(OBJS) *.dat

distclean: clean
	$(RM) -f $(EXEC)
	$(RM) *.out *.bak *~