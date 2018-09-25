#include <stdio.h>
#include <stdlib.h>
const int cols = 5;
const int rows = 5;
int **matAddition(int mat1[rows][cols], int mat2[rows][cols]);
int main(){

  int mat1[5][5] = {
    {1,1,1,2,1},
    {2,2,2,3,2},
    {3,3,3,4,3},
    {4,4,4,5,4},
    {5,5,5,4,5}};
  int mat2[5][5] = {
    {2,3,4,5,6},
    {7,8,9,10,11},
    {12,13,14,15,16},
    {17,18,19,20,21},
    {22,23,24,25,26}};

  
  int **arr = matAddition(mat1,mat2);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      printf("%d, ", arr[i][j]);
    }
    printf("\n");
  }
  
  return 0;
}

int **matAddition(int mat1[rows][cols], int mat2[rows][cols]){
  int **arr = (int **)malloc(rows*sizeof(int*));
  printf("Test1\n");
  for(int i = 0; i < rows; i++){
    arr[i] = (int *)malloc(cols*sizeof(int));
  }
  printf("Test2\n");
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      arr[i][j] = mat1[i][j]+ mat2[i][j];
    }
  }
  printf("it makes it through the function\n");
  return arr;
}
