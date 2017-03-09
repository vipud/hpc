#include <stdio.h>
#include <mpi.h>

int N;

int main(int argc, char* argv[]){
	int* a, *b, *c;
	double start_time, end_time;
	if(argc == 2){
		N = (int) argv[1];
	}
	else{
		N = 10;
	}
	a = (int *) malloc(N * sizeof(int));
	b = (int *) malloc(N * sizeof(int));
	c = (int *) malloc(N * sizeof(int));
	// init
	for(int i = 0; i < N; i++){
		a[i] = i;
		b[i] = i;
	}
	// add to c
	for(int i = 0; i < N; i++){
		c[i] = a[i] * b[i];
		#ifdef DEBUG
		printf("c[%d] = %d\n", i, c[i]);
		fflush(stdout);
		#endif
	}
	end_time = MPI_Wtime();
	printf("Time taken at size %d: %0.3f\n", N, end_time - start_time);
	free(a);
	free(b);
	free(c);
	return 0;
}
