#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "linked_list.h"
#include "../set.cvh"

#define WORD_SIZE 24
#define LO_MASK = 0x1;
#define HI_MASK = 0x008000000;
#define MASK = 0x00FFFFFF;

struct _set{
    node* head;
} 

bool add(Set* hash_set, int x){
}

bool contains(Set* hash_set, int x){
}

bool discard(Set* hash_set, int x){
}

// constructor
Set* create(int size){
    head =  create(0);
}

bool destroy(Set* hash_set){
}

// list get_all();
