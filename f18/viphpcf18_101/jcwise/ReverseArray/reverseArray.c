/* Jake Wise
 * 2018-09-25
 * Program that reverses an array. For example, {1,4,3} will become {3,4,1}. 
 * Length of array is determined by user and the values are random integers between 0 and 50.
 */

#include <stdio.h>
#include <stdlib.h>

int main () {
  int len;
  
  printf("Enter length of the array\n");
  scanf("%d", &len);

  int array[len];
  
  for (int i=0; i<len; i++) {
    array[i] = rand()%51;
  }
  
  for (int i=0; i<len; i++) {
    if (i==(len-1)) {
      printf("%d\n", array[i]);
    } else {
      printf("%d, ", array[i]);
    }
  }

  for (int i=0,j=(len-1); i<(len/2); i++, j--) {
    int temp = array[j];

    array[j] = array[i];

    array[i] = temp;

  }

  for (int i=0; i<len; i++) {
    if (i==(len-1)) {
      printf("%d\n", array[i]);
    } else {
      printf("%d, ", array[i]);
    }
  }
  
  return 0;
}
