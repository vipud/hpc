#include <stdio.h>

__global__ void hello_kernel(){
	int bid = blockIdx.x;
	int tid = threadIdx.x;
	printf("Hello from block %d, thread %d of the GPU!\n", bid, tid);
}

extern "C" void hello(){
	// do stuff here
	printf("Executing kernel...\n");
	hello_kernel<<<2,2>>>();
	cudaDeviceSynchronize();
}


