#include <stdio.h>
#include <math.h>
#include <cuda.h>

#define CHUNK_SIZE 1024
#define T unsigned long int

//make sure start is less than N/2. a is a pointer to an array of length >= N
__global__ void Fibonacci( T *a, int start) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int index = i + start;
	if (i < 2 * start - 1) {
		a[index] = (a[start - 2] * a[i]) + (a[start - 1] * a[i + 1]);
	}   
}

int main(int argc, char *argv[]) {
	int N = 0;

	if (argc == 2) {
		N = atoi(argv[1]);
	} else {
		printf("Invalid number of command line arguments.\n");
		return 1;
	}
	
	T x[3]= {1, 1, 2};
	T *d_a;

	//Allocate memory on the device
	cudaMalloc(&d_a, N*sizeof(T));
	cudaMemcpy(d_a, x, sizeof(x), cudaMemcpyHostToDevice);

	unsigned int start = 3;

	//ceiling of start - 1/1024
	while (start <= N/2 ) {
		unsigned int num_blocks = (start - 1)/CHUNK_SIZE;
		if ((start - 1) % CHUNK_SIZE != 0) {
			num_blocks++;
		}
		Fibonacci <<<num_blocks, CHUNK_SIZE>>>(d_a, start);
		start = 2 * start - 1;
	}

	//Get the results array back
	T b[N];
	cudaMemcpy(b, d_a, N*sizeof(T), cudaMemcpyDeviceToHost);

	//Print results
	for (int i = 0; i < N; i++) {
		printf("%d:\t%lu \n", i + 1, b[i]);
	}
	
	//Free device memory
	cudaFree(d_a);

	fflush(stdout);
	return 0;
}
