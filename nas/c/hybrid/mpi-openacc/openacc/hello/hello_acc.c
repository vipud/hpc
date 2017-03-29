#include <stdio.h>
#include <openacc.h>

int main(int argc, char* argv[]){
	for(int i = 0; i < 5; i++){
		printf("hello from the kernel\n!");
		fflush(stdout);
	}
}
