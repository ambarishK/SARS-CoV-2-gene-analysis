.DEFAULT_GOAL := project
CXX := g++
CXX_FLAGS := -g -std=c++17 -Werror -Wall -Wextra -march=native -ltbb

EXECUTABLES := findgenemutations.exe findrelatives.exe
CPP = $(wildcard *.cpp)
OBJ = $(CPP:%.cpp=%.o)
DEP = $(OBJ:%.o=%.d)

-include $(DEP)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@

%.exe: %.o $(OBJ)
	$(CXX) -o $@ $< $(filter-out $(EXECUTABLES:%.exe=%.o), $(OBJ)) $(CXX_FLAGS)

project: CXX_FLAGS += -O3 -DDEBUG=0
project: $(EXECUTABLES)

debug: CXX_FLAGS += -O0 -DDEBUG=1
debug: $(EXECUTABLES)

.PHONY: clean_o
clean_o:
	-rm -f $(OBJ) $(DEP)

.PHONY: clean
clean: clean_o
	-rm -f $(EXECUTABLES)
