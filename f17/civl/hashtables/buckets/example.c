#include <stdio.h>
#include <stdlib.h>

int N; // size of table
int* table;

int change(){
    t[5] = 100;
    return 1;
}

int main(int argc, char* argv[]){
    int* table; // pointer to table
    N = 10; 
    table = (int *) calloc (N, sizeof(int)); // 10 ints long 
    
    for(int i = 0; i < N; i++){
        table[i] = 10 - i;
        printf("table[%d] = %d\n", i, table[i]);
    }
    free(table);  // free memory
    return 0;
}
