#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

/* This is a hello world for omp/mpi hybrid */ 

int i = 0;

int main(int argc, char *argv[])
{
/* It's important to put this call at the begining of the program, after variable declarations. */
MPI_Init(&argc, &argv);
 
	int thread_id = 0;
	int thread_total = 1; 

 	int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int numProcs; 
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
 	
  #pragma omp parallel default(shared) private(thread_id, thread_total)
  {
  	thread_id = omp_get_thread_num();
    thread_total = omp_get_num_threads();
    printf("Hello from thread %d out of %d from process %d out of %d\n",
           thread_id, thread_total, rank, numProcs);
  }
 
	MPI_Finalize();
 
}

