#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char* argv[]){
	double start_time, end_time;
	start_time = MPI_Wtime();
	// do loop
	for(int i = 0; i < N; i++){
		// do stuff
	}
	end_time = MPI_Wtime();
	return 0;
}
