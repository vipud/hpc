#include <stdio.h>
#include <mpi.h>

int rank, size;

int main(int argc, char* argv[]){
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(rank == 0){
	printf("Hello from the CPU, running at %d cores!\n", size);
	hello();
	}
	MPI_Finalize();
	fflush(stdout);
	return 0;
}
