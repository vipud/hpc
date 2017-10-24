#include <stdio.h>
#include <stdlib.h>

int SIZE = 5;
int* table;

/* BaseHashSet... will add locks later in order to parallelize the program.
*This is just the basic adding/removing of elements to a hash set.
*Was meant to look into coarse grained hash set, but it looks as though
*coarse grained hash set just extends base hash set with a resize fxn. 
*/

// $abstract int hash_fxn(int x);
// $assume($forall (int x) hash_fxn(x) >= 0 && hash_fxn(x) < SIZE);


int hash_fxn(int x){
    return (x % 7) % SIZE;
}

int contains(int x){
	int val;
	
	val = table[hash_fxn(x)];
	return (val == x);
  //use a hash code to check and see if item is in the bucket
  
}

int add(int x){
	
	if (!contains(x)){ //if the hash set does not contain the element
		table[hash_fxn(x)] = x;
		printf("added %d to table[%d]\n", x, hash_fxn(x));
		return 1;
	}
	else{
		return 0;
	}
	 
}

int remove_element(int x){
	int val;
	
	val = table[hash_fxn(x)];
	
	if (val == x){
		int temp = val;
		table[hash_fxn(x)] = 0;
		printf("removed %d from table[%d]\n", temp, hash_fxn(x));
		return 1;
	}
	
	else{
		return 0;
	}
}

void print_hash(){
   // print table
   printf("table: ");
   for(int i = 0; i < SIZE; i++){
       printf("%d ", table[i]);
   }
   printf("\n");
}

// $input int M = 3;
// $input int INPUTS[M];

// int test0() {
//    for (int i=0; i<M; i++)
//        add(INPUTS[i]);
//    print_hash();
//    for (int i=0; i<M; i++)
//        $assert(contains(INPUTS[i]));
// }



int main(){
	table = (int *) calloc(SIZE, sizeof(int)); //table of size 5
	
	add(20);
	print_hash();
	add(17);
	print_hash();
	remove_element(17);
	print_hash();
	add(46);
	print_hash();
	add(30);
	print_hash();
	
	// test0();
	
	
	free(table);
	return 0;
}
