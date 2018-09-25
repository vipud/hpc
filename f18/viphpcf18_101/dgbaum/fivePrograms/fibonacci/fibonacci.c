#include <stdio.h>

int main(){
  printf("Enter an integer, we will print this many values of the fibonacci sequence.\n");
  int length;
  scanf("%d",&length);

  int x = 1;
  int y = 1;

  for(int i = 0; i < length; i++){
    if(i%2 == 0){
      printf("%d, ", x);
      y = y + x;
    }else{
      printf("%d, ", y);
      x = y + x;
    }
  }

  return 0;
}
