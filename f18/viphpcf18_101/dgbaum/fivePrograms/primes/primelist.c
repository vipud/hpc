#include <stdio.h>

int main(){
  int x;
  printf("Enter an integer, x, we will print the first x primes.\n");
  scanf("%d",&x);
  int a[x];
  if(x>0){
    a[0]=2;
  }else{
    printf("please enter a number > 0\n");
    return 1;
  }
  //0,1,2,3,4,5,6,7
  //2,3,5,7,11,13,17
  int counter = 3;
  for(int i = 1; i < x; i++){
    for(int j = 0; j < i; j++){
      if(counter%a[j] == 0){
	counter+=2;
	j = 0;
      }else if(j == i-1){
	a[i] = counter;
	counter+=2;
      }
    }
  }

  for(int i = 0; i < x; i++){
    printf("%d, ",a[i]);
  }
  printf("\n");
  return 0;
}
