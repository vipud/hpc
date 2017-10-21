#include <stdio.h>
#include <stdlib.h>
#include "civl.cvh"

// Concurrent Striped Hash
// Uses a fixed number of locks to write to hash table

int SIZE = 9;

struct StripedHash
{
    int* table; // table
    int* locks; // locks
};


//hash function definition
int hash(int x)
{
    return (x % SIZE);
}


// add, remove, contains
int contains(int x, struct StripedHash s)
{
    acquire(x,s);
    int slot = hash(x);
    if (s.table[slot] == x)
    {
        release(x,s);
        return 1;
    }
    else
    {
        release(x,s);
        return 0;
    }
}

int add(int x, struct StripedHash s)
{
    if(contains(x,s))
        return x;
    return 0;
}

int stripedRemove(int x, struct StripedHash s)
{
    if(contains(x,s)) //remove
        return x;
    else
        return 0;
}

void resize(struct StripedHash s)
{
    int oldCapacity = sizeof(s.table);
}

// acquire (-1), release (0) locks
void acquire(int x, struct StripedHash *s)
{
    $when(s->locks[hash(x)] > 0) *s.locks[hash(x)]--;
}

void release(int x, struct StripedHash *s)
{
    s->locks[hash(x)]++;
}


//resize policy
int policy(struct StripedHash s)
{
    return (SIZE / sizeof(s.table) > 4);
}

int main(int argc, char* argv[])
{
    struct StripedHash sHash;
    sHash.table = (int *) calloc(SIZE, sizeof(int));
    sHash.locks = (int *) calloc(SIZE, sizeof(int));
    contains(0,&sHash);
    free(sHash.table);
    free(sHash.locks);
    return 0;
}

