#include <iostream>
#include <omp.h>

int main()
{
#pragma omp parallel
    {
        int nbthreads = omp_get_num_threads();
        int rank      = omp_get_thread_num();

        std::cout << "Hello from " << rank << " among " << nbthreads << " threads." << std::endl;
    }
    return EXIT_SUCCESS;
}