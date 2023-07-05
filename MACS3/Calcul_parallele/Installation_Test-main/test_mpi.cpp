#include <iostream>
#include <mpi.h>

int main(int nargs, char *argv[])
{
    int rank, nbp;
    MPI_Init(&nargs, &argv);
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(comm, &nbp);
    MPI_Comm_rank(comm, &rank);

    std::cout << "Hello from processus " << rank << " among " << nbp << " MPI processes."
              << std::endl;

    MPI_Finalize();
    return EXIT_SUCCESS;
}