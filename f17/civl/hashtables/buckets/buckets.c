/* Aamir Majeed
 * Bucket List Hash Set
 * All items are kept in a single lock-free linked list
 * A bucket is a reference to a node in the list
 * Our array of references expands as # of buckets increase
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "buckets.h"
#include "../set.cvh"

int makeRegularKey(int x){
    int code = hash(x) & MASK;
    return reverse(code | HI_MASK);
}

int makeSentinelKey(int x){
    return reverse(x & MASK);
}

int reverse(int key){
    int loMask = LO_MASK;
    int hiMask = HI_MASK;
    int result = 0;
    for(int i = 0; i < WORD_SIZE; i++){} 
    return result;
}

bool add(Set* hash_set, int x){
}

int main(int argc, char* argv[]){

}
