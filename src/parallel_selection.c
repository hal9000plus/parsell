/*
 ============================================================================
 Name        : parallel_selection.c
 Author      : stam
 Version     :
 Copyright   : Your copyright notice
 Description : Hello MPI World in C 
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include "limits.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define ValType double
#define IS_LESS(v1, v2)  (v1 < v2)

void siftDown( ValType *a, int start, int count);

#define SWAP(r,s)  do{ValType t=r; r=s; s=t; } while(0)

void heapsort( ValType *a, int count)
{
    int start, end;

    /* heapify */
    for (start = (count-2)/2; start >=0; start--) {
        siftDown( a, start, count);
    }

    for (end=count-1; end > 0; end--) {
        SWAP(a[end],a[0]);
        siftDown(a, 0, end);
    }
}

void siftDown( ValType *a, int start, int end)
{
    int root = start;

    while ( root*2+1 < end ) {
        int child = 2*root + 1;
        if ((child + 1 < end) && IS_LESS(a[child],a[child+1])) {
            child += 1;
        }
        if (IS_LESS(a[root], a[child])) {
            SWAP( a[child], a[root] );
            root = child;
        }
        else
            return;
    }
}

int medianOfmedianSelection(int **a, int s, int e, int k){
    // if the partition length is less than or equal to 5
    // we can sort and find the kth element of it
    // this way we can find the median of n/5 partitions
    if(e-s+1 <= 5){
    	heapsort(*(a+s), e-s);
        return s+k-1;
    }
    // if array is bigger
    // we partition the array in subarrays of size 5
    // no. of partitions = n/5 = (e+1)/5
    // iterate through each partition
    // and recursively calculate the median of all of them
    // and keep putting the medians in the starting of the array
    int i=0;
    for(i=0; i<(e+1)/5; i++){
        int left = 5*i;
        int right = left + 4;
        if(right > e) right = e;
        int median = medianOfmedianSelection(a, 5*i, 5*i+4, 3);
        SWAP(*a[median], *a[i]);
    }

    // now we have array
    // a[0] = median of 1st 5 sized partition
    // a[1] = median of 2nd 5 sized partition
    // and so on till n/5
    // to find out the median of these n/5 medians
    // we need to select the n/10th element of this set (i.e. middle of it)
    return medianOfmedianSelection(a, 0, (e+1)/5, (e+1)/10);
}

int generatRandomPositiveNegitiveValue(int max , int min) {
	if ((max - 1) == RAND_MAX) {
	    return rand();
	  } else {
	    // Chop off all of the values that would cause skew...
	    long end = RAND_MAX / max; // truncate skew
	    assert (end > 0L);
	    end *= max;

	    // ... and ignore results from rand() that fall above that limit.
	    // (Worst case the loop condition should succeed 50% of the time,
	    // so we can expect to bail out of this loop pretty quickly.)
	    int r;
	    while ((r = rand()) >= end);
	    return r % max;
	  }
}
  int* generateRandomArray(int max,int min,int arraySize){
	int i=0;
	int * arr;
	arr= malloc (sizeof (int) * arraySize);
	for ( i = 0; i < arraySize; i++ )
	{
		 int temp=generatRandomPositiveNegitiveValue(max,min);
		 printf("random is %d:=%d",i,temp);
	      *(arr+i)= temp;
	}
	return arr;
}

void scatterTheArray(int sizeOfArray, int p, int arrayNumberOfPartitionElements,
		int* arr, int* receiveBuffer) {
	if (sizeOfArray % p == 0) {

		MPI_Scatter(arr, arrayNumberOfPartitionElements, MPI_INT, receiveBuffer,
				arrayNumberOfPartitionElements, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		//to implement
	}
}

void printEllementsOfProsseses(int my_rank, int arrayNumberOfPartitionElements,
		int* receiveBuffer) {
	if (my_rank != 0) {
		int ii = 0;
		printf("Process %d received data \n", my_rank);
		for (ii = 0; ii < arrayNumberOfPartitionElements; ii++) {
			printf("element %d=%d process %d ", ii, (receiveBuffer[ii]),
					my_rank);
			printf("\n");
		}
	} else {
		int ii = 0;
		printf("Process %d received data \n", my_rank);
		for (ii = 0; ii < arrayNumberOfPartitionElements; ii++) {
			printf("element %d=%d process %d ", ii, (receiveBuffer[ii]),
					my_rank);
			printf("\n");
		}
	}
}

  /*
   * input arguments argv[1]=size of random array argv[2]=the rank of the ellement we are looking for
   */
int main(int argc, char* argv[]){
	int my_rank; /* rank of process */
	int p;       /* number of processes */
	int source;   /* rank of sender */
	int dest;     /* rank of receiver */
	int tag=0;    /* tag for messages */
	int sizeOfArray,numberOfEllementsLeft,constantParametrized=100;
	MPI_Status status ;   /* return status for receive */
	int *sendBuffer,*receiveBuffer,*arr;
	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//initialiate after the call to MPI_INIT to get the correct argv ....
	sizeOfArray=atoi(argv[1]);
	int arrayNumberOfPartitionElements=sizeOfArray/p;
	receiveBuffer=(int *)malloc(arrayNumberOfPartitionElements*sizeof(int));/// all allocate a place holder for scattered data
	//int kSmallestEllement=atoi(argv[2]);
	if (my_rank !=0){//slaves ---usually black :0

	}
	else{//master
		arr=(int *)malloc(sizeOfArray*sizeof(int));//dont waste memmory ..only master allocates the large array
		arr=generateRandomArray(100,INT_MIN,sizeOfArray);
		numberOfEllementsLeft=sizeOfArray;//INITIAL VALUE OF REMAINIGN NUMBER OF ELLEMENTS
	}
	scatterTheArray(sizeOfArray, p, arrayNumberOfPartitionElements, arr,
			receiveBuffer);
	printEllementsOfProsseses(my_rank, arrayNumberOfPartitionElements,
			receiveBuffer);//tracing

	while( numberOfEllementsLeft < (constantParametrized*p) ){
		int n=(sizeof(receiveBuffer)/sizeof(receiveBuffer[0]));
		int mom=medianOfmedianSelection(&receiveBuffer,n-1,0,n/2);
		printf("THE MON ISISISISIS :+%d",mom);
		break;
	}

	MPI_Finalize(); 
	
	
	return 0;
}





