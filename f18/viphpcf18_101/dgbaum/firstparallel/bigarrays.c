#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>

int main(){
  int size = 10000; //this is the largest I can make the arrays...
                  //otherwise I get segmentation fault (core dump)
  int **array1 = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    array1[i] = (int *)malloc(size * sizeof(int));


  int **array2 = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    array2[i] = (int *)malloc(size * sizeof(int));

  int **sumArray = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    sumArray[i] = (int *)malloc(size * sizeof(int));
  
  int **sumArray2 = (int **)malloc(size * sizeof(int *));
  for(int i = 0; i < size; i++)
    sumArray2[i] = (int *)malloc(size * sizeof(int));
  

  srand((int) time(NULL));

  
  
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      array1[i][j] = rand()%10000;
      array2[i][j] = rand()%10000;
    }
  }

  //here we will do the normal adding
  double start,end;
  double time_elapsed;

  start = omp_get_wtime();

  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      sumArray[i][j] = array1[i][j] + array2[i][j];
    }
  }

  end = omp_get_wtime();

  time_elapsed = ((double) (end - start));

  printf("Normal adding took %f seconds.\n",time_elapsed);
  
  start = omp_get_wtime();
  //here we will do the parallel adding
  int i,j;
#pragma omp parallel for default(none) shared(sumArray2,size,array1,array2) private(i,j) collapse(2)
  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      sumArray2[i][j] = array1[i][j] + array2[i][j];
    }
  }

  end = omp_get_wtime();

  time_elapsed = ((double) (end - start));

  printf("Adding with parallelism took %f seconds.\n",time_elapsed);

  //verifying that the parallel adding is adding properly.
  bool same = true;
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      if(sumArray[i][j] != sumArray2[i][j]){
	same = false;
	i = j = size;
	break;
      }
    }
  }
  if(same){
    printf("Our parallel code is producing the proper solution.\n");
  }else{
    printf("Our parallel code does something wrong and doesn't work properly.\n");
  }
  return 0;
}
