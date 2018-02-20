#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]){
    MPI_Init(NULL, NULL); // initialize MPI
    int n_procs; // number of processors
    int n_rank;  // current processor rank
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs); // get number of processors, store in n_procs
    MPI_Comm_rank(MPI_COMM_WORLD, &n_rank); // get current processor, store in n_rank
    // print from each processor
    printf("Hello world from processor %d of %d.\n", n_rank, n_procs);
    MPI_Finalize();
    return 0;
}
