#include <stdio.h>
#include <stdlib.h>

int N; // size of table
int* table; // pointer to table

int add(int x){
    return 1;
}

int contains(int x){
    return 1;
}

int remove(int x){
    return 1;
}

int main(int argc, char* argv[]){
    N = 10; 
    table = (int *) malloc (N * sizeof(int)); // 10 ints long 
    // do stuff
     
    free(table);  // free memory
    return 0;
}
