#include <stdio.h>
#include <stdlib.h>

// Concurrent Striped Hash - Single core
// Uses a fixed number of locks to write to hash table

int LIMIT = 9;

struct StripedHash{
    int* table; // table
    int* locks; // locks
};

int hash(int x){
    return (x % 11) % 9;
}

// add,remove,contains
int contains(int x, struct StripedHash s){
    return x;
}

int add(int x, struct StripedHash s){
    if(contains(x,s))
        return x;
    return 0;
}

int remove(int x, struct StripedHash s){
    if(contains(x,s)) //remove
        return x;
    else
        return 0;
}

// acquire, release locks

int acquire(int x){
    return x;
}

int release(int x){
    return x;
}

int main(int argc, char* argv[]){
    struct StripedHash s;
    return 0;
}

