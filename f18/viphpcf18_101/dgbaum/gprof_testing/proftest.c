#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int size = 10000;
void serialFillArray(int **arr);
void parallelFillArray(int **arr);
void serialSumArray(int **sum,int **arr1, int **arr2);
void parallelSumArray(int **arr);
int main(){
  srand((int)time(NULL));
  int **array1 = (int **)malloc(size*sizeof(int *));
  int **array2 = (int **)malloc(size*sizeof(int *));
  int **sumArray1 = (int **)malloc(size*sizeof(int *));
  int **sumArray2 = (int **)malloc(size*sizeof(int *));
  for(int i = 0; i < size;i++){
    array1[i]=(int *)malloc(size*sizeof(int));
    array2[i]=(int *)malloc(size*sizeof(int));
    sumArray1[i]=(int *)malloc(size*sizeof(int));
    sumArray2[i]=(int *)malloc(size*sizeof(int));
  }
  printf("%d\n",size);

  serialFillArray(array1);
  serialFillArray(array2);
  serialSumArray(sumArray1,array1,array2);
  return 0;
}
void serialFillArray(int **arr){
  for(int i=0;i<size;i++){
    for(int j = 0;j<size;j++){
      arr[i][j] = rand()%10;
    }
  }
}
void parallelFillArray(int **arr){

}
void serialSumArray(int **sum, int **arr1, int **arr2){

  for(int i = 0; i < size;i++){
    for(int j = 0; j < size;j++){
      sum[i][j] = arr1[i][j] + arr2[i][j];
    }
  }
  
}
void parallelSumArray(int **arr){

}
