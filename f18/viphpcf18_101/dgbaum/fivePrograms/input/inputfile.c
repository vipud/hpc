#include <stdio.h>

int main(int argc, char *argv[]){
  FILE *fp;
  char *filename;
  char ch;

  if(argc < 2){
    printf("Missing Filename\n");
    return 1;
  }else{
    filename = argv[1];
    printf("Filename : %s\n", filename);
  }

  fp = fopen(filename,"r");


  if(fp){
    printf("File Contents:\n");

    while((ch = fgetc(fp)) != EOF){
      printf("%c",ch);
    }

  }else{
    printf("We failed to open the file\n");   
  }
  return 0;
}
