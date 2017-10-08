#include <stdio.h>
#include <stdlib.h>

// Cuckoo Hash

int LIMIT = 9;

struct CuckooHash{
    int* table0; // uses hash 0
    int* table1; // uses hash 1
};

int hash0(int x){
    return (x % 11) % 9;
}

int hash1(int x){
    // printf("%d % 13 % 9 = %d\n", x, (x % 13) % 9); 
    return (x % 13) % 9;
}


int swap_hash(int* x1, int* x2){
    //printf("swap_hashping %d and %d\n", *x1, *x2);
    int temp = *x2;
    //printf("temp = %d\n", temp);
    *x2 = *x1;
    *x1 = temp;
    //printf("x1 is %d, x2 is %d\n", *x1, *x2);
    return *x2;
}


// check if element exists in either table
int contains(int x, struct CuckooHash h){
    //printf("hash0(%d) = %d\n", x, hash0(x));
    //printf("hash1(%d) = %d\n", x, hash1(x));
    int val_t0, val_t1;
    
    //printf("Table 0 [%d] = %d\n", hash0(x), h.table0[hash0(x)]);
    val_t0 = (h.table0[hash0(x)] == x);
    //printf("Table 1 [%d] = %d\n", hash1(x), h.table1[hash1(x)]);
    //printf("val_t0 = %d\n", val_t0);
    val_t1 = (h.table1[hash1(x)] == x);
    //printf("val_t1 = %d\n", val_t1);
    return val_t0 || val_t1;

}

// add a value to the hash
int add(int x, struct CuckooHash h){
    if(contains(x, h)){
        printf("%d exists in the set.\n", x);
        return 0;
    }
    for(int i = 0; i < LIMIT; i++){
        if(swap_hash(&(h.table0[hash0(x)]), &x) == 0){
            //printf("swap_hashped 0\n");
            return 1;
        } 
        //printf("swap_hashping %d with %d\n", h.table1[hash1(x)], x);
        if (swap_hash(&(h.table1[hash1(x)]), &x) == 0){
            return 1;
        }
    }
   // printf("did not add %d", x);
}

// print currently stored values
void print_hash(struct CuckooHash h){
    // print table0
    printf("table0: ");
    for(int i = 0; i < LIMIT; i++){
        printf("%d ", h.table0[i]);
    }
    // print table1
    printf("table1: ");
    for(int i = 0; i < LIMIT; i++){
        printf("%d ", h.table1[i]);
    }
    printf("\n\n");
}


int main(){
    struct CuckooHash h;

    // initialize
    h.table0 = (int *) calloc(9, sizeof(int));
    h.table1 = (int *) calloc(9, sizeof(int));

    print_hash(h);
    printf("adding 1\n");
    add(1, h); 
    print_hash(h);
    printf("adding 2\n");
    add(3,h);
    printf("adding 134\n");
    add(135,h);
    printf("adding 133\n");
    add(133, h);
    printf("adding 66\n");
    add(66, h);
    print_hash(h);
    printf("adding 45\n");
    add(45, h);
    print_hash(h);
    add(100,h);
    print_hash(h);
    add(254,h);
    print_hash(h);
    add(101,h);
    print_hash(h);
    add(2,h);
    print_hash(h);
    add(103,h);
    add(104,h);
    add(105,h);
    add(106,h);
    add(107,h);
    add(108,h);
    add(109,h);
    add(110,h);
    print_hash(h);
    add(111,h);
    print_hash(h);
    add(112,h);
    print_hash(h);
    add(113,h);
    print_hash(h);
    add(134, h);
    print_hash(h);
    if(contains(44, h)) printf("44 exists\n");
    if(contains(45, h)) printf("45 exists\n"); 
    if(contains(1, h)) printf("1 exists\n"); 
    print_hash(h);

    free(h.table0);
    free(h.table1);
    return 0;

}
