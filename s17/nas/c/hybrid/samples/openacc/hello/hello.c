#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]){
	double start, end;
	start = MPI_Wtime();
	printf("Hello world!\n");
	end = MPI_Wtime();
	printf("Time taken, %0.3f\n",end - start);
	return 0;
}
