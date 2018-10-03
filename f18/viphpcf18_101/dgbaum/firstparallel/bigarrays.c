#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int main(){
  int size = 10000; //this is the largest I can make the arrays...
                  //otherwise I get segmentation fault (core dump)

  int **arr = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    arr[i] = (int *)malloc(size * sizeof(int));
  
  
  int **array1 = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    array1[i] = (int *)malloc(size * sizeof(int));


  int **array2 = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    array2[i] = (int *)malloc(size * sizeof(int));

  int **sumArray = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    sumArray[i] = (int *)malloc(size * sizeof(int));
  
  
  

  srand((int) time(NULL));

  
  
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      array1[i][j] = rand()%10000;
      array2[i][j] = rand()%10000;
    }
  }

  //here we will do the normal adding
  clock_t start,end;
  double cpu_time_used;

   start = clock();

  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      sumArray[i][j] = array1[i][j] * array2[i][j];
    }
  }

   end = clock();

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("Normal adding took %f seconds.\n",cpu_time_used);
  
  start = clock();
  //here we will do the parallel adding
  int i,j;
#pragma omp parallel for default(none) shared(sumArray,size,array1,array2) private(i,j)
  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      sumArray[i][j] = array1[i][j] * array2[i][j];
    }
  }

  end = clock();

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("Adding with parallelism took %f seconds.\n",cpu_time_used);
  return 0;
}
