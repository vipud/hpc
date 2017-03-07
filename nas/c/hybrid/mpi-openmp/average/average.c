#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MASTER 0
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

double sum(double* start, int len);
int main(int argc, char *argv[]){

	int numtasks;
	int taskid;
	int rc;
	double totalSum;

	double arr[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
	double* parr = arr;

	int arrSize = NELEMS(arr);

	//double *meme = arr+2;
	//printf("%f", *meme);

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	//printf("MPI task %d has started... \n", taskid);

	int splitLength = arrSize/numtasks;
	int leftover = arrSize%numtasks;
	double thisSum;

	if(taskid == (numtasks-1)){
		thisSum = sum(parr+(splitLength*taskid), splitLength+leftover);
		printf("Sum for thread %d:  %f\n",taskid,thisSum);
	}else{
		thisSum = sum(parr+(splitLength*taskid),splitLength);
	        printf("Sum for thread %d:  %f\n",taskid,thisSum);
	}
	//double thisAvg = avg(parr+(splitLength*taskid),splitLength);

	//double thisAvg = avg(parr+(2*taskid),2);
	//printf("%f\n",thisAvg);

	rc = MPI_Reduce(&thisSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	if(rc != MPI_SUCCESS)
		printf("%d: falure on mpc_reduce\n",taskid);

	//printf("Average is equal to %f \n", avg(parr+(2*taskid),2));


	if(taskid == MASTER){
		//printf("Average Sum %f", averageSum);
		double finalAverage = totalSum/arrSize;
		printf("Average sum: %f\n",totalSum);
		printf("The final average is %f\n", finalAverage);
	}

	MPI_Finalize();
	return 0;
}

double sum(double* start, int len){
	double total = 0;
	int i;
	for(i = 0; i < len; i++){
		total+=*(start+i);
	}
	return total;
	//return (*start+*(start+1))/2;
}
