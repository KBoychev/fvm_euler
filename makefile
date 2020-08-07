

COMPILER = g++
COMPILER_FLAGS = -O2 -std=c++11 -flto -fopenmp 

SOURCES = $(wildcard ./includes/*.cpp)
OBJECTS = $(patsubst ./includes/%.cpp,./bin/%.o,$(SOURCES))
 
all: euler_solver clean
	
euler_solver: $(OBJECTS) ./bin/solver.o
	$(COMPILER) $(COMPILER_FLAGS) $(OBJECTS) ./bin/solver.o -o ./bin/euler_solver
./bin/%.o:./includes/%.cpp
	$(COMPILER) $(COMPILER_FLAGS) -I./includes -c $< -o $@
./bin/solver.o:./solver/solver.cpp
	$(COMPILER) $(COMPILER_FLAGS) -I./includes -c ./solver/solver.cpp -o ./bin/solver.o
clean:
	rm -f ./bin/*.o
