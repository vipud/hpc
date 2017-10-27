/* Aamir Majeed
 * Bucket List Hash Set
 * All items are kept in a single lock-free linked list
 * A bucket is a reference to a node in the list
 * Our array of references expands as # of buckets increase
 */

#include <stdio.h>
#include <stdlib.h>
#include "linkedlist.h"

int N; // size of table, always power of 2

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
