CC = mpicc
CXX = mpic++
LIB = -lmpi
MPIOPT = -oversubscribe
#include Make_msys2.inc
#include Make_osx.inc

CXXOPT = -g -Wall -pedantic -std=c++11 -fopenmp

all:	test_mpi.exe test_omp.exe test

test_mpi.exe:	test_mpi.cpp
	$(CC) $(CXXOPT) -o test_mpi.exe test_mpi.cpp $(LIB)

test_omp.exe:	test_omp.cpp
	$(CC) $(CXXOPT) -o test_omp.exe test_omp.cpp $(LIB)

test: test_mpi.exe test_omp.exe
	mpiexec $(MPIOPT) -np 2 ./test_mpi.exe
	mpiexec $(MPIOPT) -np 16 ./test_mpi.exe
	./test_omp.exe

clean:
	@rm -f *.exe *~ 
