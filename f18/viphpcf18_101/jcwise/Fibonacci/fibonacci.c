/* Jake Wise
 * 2018-09-25
 * Program that computes Fibonacci sequence, starting at 0, up to the given number.
 */

#include <stdio.h>

int main() {
  int input;
  long num1 = 0;
  long num2 = 1;
  long currentFib;
  long temp;

  printf("Enter the highest Fibonacci number you want.\n");
  scanf("%d", &input);
  printf("The Fibonacci series up to %d is:\n", input);
  
  if(input==0) {
    printf("0\n");
    return 0;
  }
  if(input==1) {
    printf("0, 1\n");
    return 0;
  }
  printf("0, 1, ");  
  for (int i=0; i < input-1; i++) {
    currentFib = num1 + num2;
    temp = num2;
    num2 = currentFib;
    num1 = temp;
    
    if (i==(input-2)) {
      printf("%ld\n", currentFib);
    } else {
    printf("%ld, ", currentFib);
    }
  }
  return 0;
}
