#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "linked_list.h"
#include "../set.cvh"

#define WORD_SIZE 24
#define LO_MASK  0x1
#define HI_MASK 0x008000000
#define MASK  0x00FFFFFF

int size;
struct _window{
    node* pred;
    node* curr;
};
typedef struct _window Window;

Window* find(node* head, int key){
    node* pred = head;
    node* curr = head->next;
    Window* new_window;
    while(curr->val < key){
        pred = curr;
        curr = pred->next;
    }
    new_window->pred = pred;
    new_window->curr = curr;
    return new_window;
};

int hash(int key){
    return key % size;
};

int hashCode(int key){
    return hash(key) & MASK; // restricted by mask 
}

int reverse(int key);

int makeRegularKey(int key);

int makeSentinelKey(int key);

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

/* To insert a key into our hash table
 * The key must be hashed using recursive split-ordering
 * ...following the pointer to the appropriate location in the sorted items list
 * ...and traversing the list until the key's proper location is found 
 */ 


