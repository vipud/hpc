
#include <time.h>
#include <stdio.h>

int main() {
  int row;
  int col;

  printf("Enter the number of rows: ");
  scanf("%d", &row);
  printf("Enter the number of columns: ");
  scanf("%d", &col);

  int mat1[row][col];
  int mat2[row][col];
  int matSum[row][col];
  double start, finish;
  printf("Enter the values for matrix 1 (press enter after each number)\n");
 
  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      scanf("%d", &mat1[i][j]);
    }
  }

  printf("Enter the values for matrix 2\n");

  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      scanf("%d", &mat2[i][j]);
    }
  }
  start = omp_get_wtime(); 
#pragma omp parallel for
  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      matSum[i][j] = mat1[i][j] + mat2[i][j];
    }
  }
  finish = omp_get_wtime();

  double time = finish - start;

  printf("First matrix:\n");
  printf("time: %d\n", time);
  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      printf("%d ", mat1[i][j]);
    }
    printf("\n");
  }

  printf("Second matrix:\n");
 
  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      printf("%d ", mat2[i][j]);
    }
    printf("\n");
  }

  printf("Sum matrix:\n");
  
  for(int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      printf("%d ", matSum[i][j]);
    }
    printf("\n");
  }

  return 0;
}
