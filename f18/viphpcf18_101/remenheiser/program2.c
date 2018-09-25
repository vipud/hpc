// check if a number is prime or not

#include <stdio.h>

int isPrime(int x);

int main(){
  int input;

  printf("Enter a positive integer: \n");
  scanf("%d", &input);
 
  if (isPrime(input)){
      printf("%d is prime \n", input);
    }
  else {
    printf("%d is not prime \n", input);
  }
  return 0;
}

int isPrime(int x){
  
  if (x == 2){
    return 1;
  }
  
  for (int i = 2; i < x/2; i++){
    if (x % i == 0){
      return 0;
    }
  }  
  return 1;
}
