#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "linked_list.h"
#include "../set.cvh"

#define WORD_SIZE 24
#define LO_MASK = 0x1;
#define HI_MASK = 0x008000000;
#define MASK = 0x00FFFFFF;
int size;

/* A bucket consists of a reference into our list
 * When resizing, move buckets rather than individual items
 * Implemented using: 
 * - The combination of a modulo-size hash and a 2^i table size.
 * - Recursive split-ordering
 *
 * Given i, where 2^i = SIZE, and b, a logical table bucket
 * An item k maintains 'k mod 2^i = b' 
 * When resizing, such that SIZE = 2^(i+1),
 * the items in the bucket are split between bucket b... 
 * ...and a new bucket, such that'k mod 2^(i+1) = b + 2^i'.
 * Resize is incremental
 *
 * Recursive split-ordering is achieved through binary reversal
 * - reversing most-significant bits with least-significant bits
 * 00011 -> 11000
 *
 * Used in (concurrent) LockFreeHashSet, holds an array of bucketlists
 */

/* Atomic Markable Reference
 * compareAndSet(expected-val,new-val,expected-mark,new-mark) 
 * updates if curr=expected
 * in Shavit's implementation, checks node vals
 * add: if nodes not equal, stop and don't add
 * otherwise continue
 */

struct _set{
    node* head;
}  

struct _window{
    node* pred;
    node* curr;
}

typedef struct _set Set;
typedef struct _window Window;

int reverse(int key);

int makeKey(int key);

int makeSentinelKey(int key);


/* To insert a key into our hash table
 * The key must be hashed using recursive split-ordering
 * ...following the pointer to the appropriate location in the sorted items list
 * ...and traversing the list until the key's proper location is found 
 */ 
bool add(Set* hash_set, int x);

bool contains(Set* hash_set, int x);

bool discard(Set* hash_set, int x);

Set* create(int size);

bool destroy(Set* hash_set);
