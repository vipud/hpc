// Cuckoo Hash

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int LEN = 5; // length of each table
struct CuckooHash h;

struct CuckooHash{
   int* table0; // uses hash 0
   int* table1; // uses hash 1
   int items; // number of items
};
//int hash0(int x){
//    return (x % 11) % LEN;
//}
//
//int hash1(int x){
//    return (x % 13) % LEN;
//}

$abstract int hash0(int x);
$assume($forall (int x) hash0(x) >= 0 && hash0(x) < LEN);
$abstract int hash1(int x);
$assume($forall (int x) hash1(x) >= 0 && hash1(x) < LEN);

$assume($forall (int x)
        ($forall (int y)
          ((x == y) => hash0(x) != hash1(y))
        )
       );
       
//Bounding to eliminate error created by abstract functions
$assume($forall (int x)
		 ($forall (int y)
		  ($forall (int z)
		  	((hash0(x) == hash0(y) && hash0(x) == hash0(z)) =>
		  		!(hash1(x) == hash1(y) && hash1(x) == hash1(z)))
		  )
		 )
		);

int swap_hash(int* x1, int* x2){
   int temp = *x2;
   *x2 = *x1;
   *x1 = temp;
   return *x2;
}

// check if element exists in either table
int contains(int x, struct CuckooHash h){
   int val_t0, val_t1;
   
   val_t0 = (h.table0[hash0(x)] - x);
   val_t1 = (h.table1[hash1(x)] - x);
   return (val_t0 * val_t1) == 0;
   
}

// checks # items / both table lengths
int policy(struct CuckooHash h){
    double load_factor = h.items / (2*LEN);
    return load_factor > 0.75;
}

void resize(struct CuckooHash h){
    // store previous hash values
    int* temp_table0 = (int*)malloc(LEN * sizeof(int));
    int* temp_table1 = (int*)malloc(LEN * sizeof(int));
    memcpy(temp_table0, h.table0, sizeof(h.table0));    
    memcpy(temp_table1, h.table1, sizeof(h.table0));    
    // double size in table0, table1
    LEN *= 2;
    free(h.table0);
    free(h.table1);
    // empty tables
    h.table0 = (int *) calloc(LEN, sizeof(int)); 
    h.table1 = (int *) calloc(LEN, sizeof(int));
    // rehash 
}

// add a value to the hash
int add(int x, struct CuckooHash h){
   if(contains(x, h)){
       printf("%d exists in the set.\n", x);
       return 0;
   }
   for(int i = 0; i < LEN; i++){
       if(swap_hash(&(h.table0[hash0(x)]), &x) == 0)
           h.items += 1;
           return 1;
       if (swap_hash(&(h.table1[hash1(x)]), &x) == 0)
           h.items += 1;
           return 1;
   }
   resize();
   add(x,h);
   return;
}

int remove(int x, struct CuckooHash h){
   int val_t0, val_t1;

   val_t0 = h.table0[hash0(x)];
   val_t1 = h.table1[hash1(x)];
   if(val_t0 == x){
        int temp = val_t0;
        h.table0[hash0(x)] = 0;
        printf("removed %d from table0[%d]\n", temp, hash0(x));
        return 1;
   }
   else if(val_t1 == x){
        int temp = val_t1;
        h.table0[hash0(x)] = 0;
        printf("removed %d from table0[%d]\n", temp, hash0(x));
        return 1;
   }
   else{
       return 0;
   }
}

// print currently stored values
void print_hash(struct CuckooHash h){
   // print table0
   printf("table0: ");
   for(int i = 0; i < LEN; i++){
       printf("%d ", h.table0[i]);
   }
   printf("\n");
   // print table1
   printf("table1: ");
   for(int i = 0; i < LEN; i++){
       printf("%d ", h.table1[i]);
   }
   printf("\nLoad factor = %d/%d = %d", h.items, LEN, (double)h.items/LEN);
}

int test1(){
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

$input int M = 5;
$input int INPUTS[M];
$assume($forall (int x : 0 .. M-1) INPUTS[x] > 0);
$assume($forall (int x : 0 .. M-1)
        ($forall (int y : 0 .. M-1)
         (x != y) => (INPUTS[x] != INPUTS[y])
         ));

int test2() {
   for (int i=0; i<M; i++)
       add(INPUTS[i], h);
   print_hash(h);
   for (int i=0; i<M; i++)
       $assert(contains(INPUTS[i], h));
}
int main() {
   // initialize
   h.table0 = (int *) calloc(LEN, sizeof(int));
   h.table1 = (int *) calloc(LEN, sizeof(int));
   
   // test1();
   test2();
   
   free(h.table0);
   free(h.table1);
}
