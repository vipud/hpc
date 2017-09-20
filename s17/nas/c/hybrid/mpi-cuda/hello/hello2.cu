#include <stdio.h>
#include <mpi.h>
#include <cuda.h>

__global__ void hello(){
	printf("Hello from the kernel!");
}

int main(int argc, char* argv[]){
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("Hello from the processor!");
	hello<<<1,1>>>();
	MPI_Finalize();
	return 0;
}

