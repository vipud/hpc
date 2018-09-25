#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

  int nthreads, tid;

  /*get a team of threads and give them their own copies of variables*/
#pragma omp parallel private(nthreads, tid)
  {
  tid = omp_get_thread_num();
  printf("Hello world from thread: %d\n", tid);

  if (tid == 0){
    nthreads = omp_get_num_threads();
    printf("Number of threads - %d\n",nthreads);
  }
  
  }/*All threads join master thread and disband */

  
}
