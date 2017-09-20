#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <openacc.h>

#define N 10000

int main(int argc, char* argv[]){
	double start_time, end_time;
	start_time = MPI_Wtime();
	// do loop
    int x = 0;
    #pragma acc kernels
	for(int i = 0; i < N; i++){
		// do stuff
        x += i * i;
	}
	end_time = MPI_Wtime();
    printf("Time taken = %0.2f, x = %d.\n", (end_time - start_time), x);
	return 0;
}
