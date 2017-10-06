#include <stdio.h>
#include <stdlib.h>

// Cuckoo Hash

int LIMIT = 9;

struct CuckooHash{
    int** table0; // uses hash 0
    int** table1; // uses hash 1
};

int hash0(int x){
    return (x % 11) % 9;
}

int hash1(int x){
    return (x % 13) % 9;
}

int swap(int* x1, int* x2){
    printf("swapping %d and %d\n", *x1, *x2);
    int temp = *x1;
    *x2 = *x1;
    *x1 = temp;
    printf("x1 is %d, x2 is %d\n", *x1, *x2);
}

// check if element exists in either table
int contains(int x, struct CuckooHash h){
    printf("hash0(%d) = %d\n", x, hash0(x));
    printf("hash1(%d) = %d\n", x, hash1(x));
    int val_t0, val_t1;
    
    val_t0 = (*(h.table0[hash0(x)]) == x);
    printf("val_t0 = %d\n", val_t0);
    val_t1 = (*(h.table1[hash1(x)]) == x);
    printf("val_t1 = %d\n", val_t1);
    return val_t0 || val_t1;

}

// add a value to the hash
int add(int x, struct CuckooHash h){
    if(contains(x, h)){
        printf("%d exists in the set.\n", x);
        return 0;
    }
    for(int i = 0; i < LIMIT; i++){
        printf("at %d\n", i); 
        if((x == swap(h.table0[hash0(x)], &x)) == 0){
            printf("swapped 0\n");
            return 1;
        } else if ((x == swap(h.table1[hash1(x)], &x) == 0)){
            return 1;
        }
    }
    add(x,h);
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
    h.table0 = (int **) calloc(9, sizeof(int*));
    h.table1 = (int **) calloc(9, sizeof(int*));

    print_hash(h);
    add(1, h); 
    print_hash(h);
    add(67, h);
    print_hash(h);
    add(45, h);
    print_hash(h);
    if(contains(44, h)) printf("44 exists\n");
    if(contains(45, h)) printf("45 exists\n"); 
    print_hash(h);

    free(h.table0);
    free(h.table1);
    return 0;

}
