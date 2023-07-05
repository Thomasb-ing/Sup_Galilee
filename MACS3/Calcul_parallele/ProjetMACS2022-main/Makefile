# include à modifier selon sa configuration
CXX = mpic++
CXXFLAGS = -std=c++17 -fopenmp
ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG -fopenmp
else
CXXFLAGS2 = ${CXXFLAGS} -O2 -march=native -Wall -fopenmp
CXXFLAGS += -O3 -march=native -Wall
endif
LIB=-lsfml-graphics -lsfml-window -lsfml-system


ALL= vortexSimulation.exe

default:	help
all: $(ALL)

clean:
	@rm -fr *.o *.exe *~ *.png

OBJS= vortex.o screen.o runge_kutta.o cloud_of_points.o cartesian_grid_of_speed.o vortexSimulation.o

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ $^

vortexSimulation.exe: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIB)

help:
	@echo "Available targets : "
	@echo "    all                           : compile all executables"
	@echo "    syracuse_simple.exe           : compile simple syracuse sequence executable"
	@echo "    syracuse_orbite.exe           : compile syracuse with orbit executable"
	@echo "    simple_automate_cellulaire.exe: compile simple cellular automata executable"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
