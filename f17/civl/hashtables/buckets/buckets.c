/* Aamir Majeed
 * Bucket List Hash Set
 * All items are kept in a single lock-free linked list
 * A bucket is a reference to a node in the list
 * Our array of references expands as # of buckets increase
 */

#include <stdio.h>
#include <stdlib.h>
#include "linked_list.h"
#include "../set.cvh"

Set* create(int size){
    Set* s;
    hs = 
    hs->table = (int *) malloc(size * sizeof(int));
    return hs;
}

int N; // size of table, always power of 2

int add(int x){
    return 1;
}

int contains(int x){
    return 1;
}

int n_remove(int x){
    return 1;
}

int main(int argc, char* argv[]){
    N = 10; 
    int* table;
    bucket_list* bucket; // malloc
    bucket = (bucket_list*) malloc(sizeof(bucket_list));
    bucket->val = 10;
    table = (int *) malloc (N * sizeof(int)); // 10 ints long 
    // do stuff
    printf("bucket->val = %d\n", bucket->val); 
    free(table);  // free memory
    return 0;
}
