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
#include <stddef.h>
#include <time.h>
#define SWAPINT(r,s)  do{int t=r; r=s; s=t; } while(0)
#define SWAPDOUBLE(r,s)  do{double t=r; r=s; s=t; } while(0)
#define NEWARRAY(x, n) do { (x) = calloc((n), sizeof *(x)); } while (0)
#define M 7
#define NSTACK 50
#define NR_END 1
#define FREE_ARG char*

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _b : _a; })

typedef struct LEG{
	int *L;
	int *E;
	int *G;
} LEG;

typedef struct LOCALCOUNTS{
	int li;
	int ei;
	int gi;
} LOCALCOUNTS;


LEG* New_LEG(int size){
	 LEG *newLEG = malloc(sizeof(LEG));
	    assert(newLEG != NULL);
	NEWARRAY(newLEG->L,size);
	NEWARRAY(newLEG->E,size);
	NEWARRAY(newLEG->G,size);
	return newLEG;
}



void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

void inline quicksort(unsigned long n, int arr[]){
	if(n>1){
	unsigned long i,ir=n,j,k,l=1,*istack;
	int jstack=0;
	int a,temp;
	istack=lvector(1,NSTACK);
	for (;;) {
		//Insertion sort when subarray small enough.
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			//Pop stack and begin a new round of parti-
			l=istack[jstack--];
			//tioning.
		} else {
			k=(l+ir) >> 1;
			//Choose median of left, center, and right el-
			SWAPINT(arr[k],arr[l+1]);
			//ements as partitioning element a. Also
			if (arr[l] > arr[ir]) {
				//rearrange so that a[l] ≤ a[l+1] ≤ a[ir].
				SWAPINT(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				SWAPINT(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				SWAPINT(arr[l],arr[l+1]);
			}
			i=l+1;
			//Initialize pointers for partitioning.
			j=ir;
			a=arr[l+1];
			//Partitioning element.
			for (;;) {
				//Beginning of innermost loop.
				do i++; while (arr[i] < a);
				//Scan up to find element > a.
				do j--; while (arr[j] > a);
				//Scan down to find element < a.
				if (j < i) break;
				//Pointers crossed. Partitioning complete.
				SWAPINT(arr[i],arr[j]);
				//Exchange elements.
			}
			//End of innermost loop.
			arr[l+1]=arr[j];
			//Insert partitioning element.
			arr[j]=a;
			jstack += 2;
			//Push pointers to larger subarray on stack, process smaller subarray immediately.
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_lvector(istack,1,NSTACK);
	}else{
			printf("SORT ARRAY trying to short an array with n=%d",n);
			exit(1);
		}
}

/*Returns the k th smallest value in the array arr[1..n] . The input array will be rearranged
to have this value in location arr[k] , with all smaller elements moved to arr[1..k-1] (in
		arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).*/
int inline selectkthelem(unsigned long k, unsigned long n, int *arr){
	if(n>1 && k<=n  ){
	unsigned long i,ir,j,l,mid;
	int a,temp;
	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			//Active partition contains 1 or 2 elements.
			if (ir == l+1 && arr[ir] < arr[l]) {
				//Case of 2 elements.
				SWAPINT(arr[l],arr[ir]);
			}
			if(k<=0){
				printf("the index is -1 ");
				exit(1);
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			//Choose median of left, center, and right el-
			SWAPINT(arr[mid],arr[l+1]);
			//ements as partitioning element a. Also
			if (arr[l] > arr[ir]) {
				//rearrange so that arr[l] ≤ arr[l+1],
				SWAPINT(arr[l],arr[ir]);
				//arr[ir] >= arr[l+1];
			}
			if (arr[l+1] > arr[ir]) {
				SWAPINT(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				SWAPINT(arr[l],arr[l+1]);
			}
			i=l+1;
			///Initialize pointers for partitioning.
			j=ir;
			a=arr[l+1];
			//Partitioning element.
			for (;;) {
				//Beginning of innermost loop.
				do i++; while (arr[i] < a);
				//Scan up to find element > a.
				do j--; while (arr[j] > a);
				//Scan down to find element < a.
				if (j < i) break;
				//Pointers crossed. Partitioning complete.
				SWAPINT(arr[i],arr[j]);
			}
			//End of innermost loop.
			arr[l+1]=arr[j];
			//Insert partitioning element.
			arr[j]=a;
			if (j >= k) ir=j-1;
			//Keep active the partition that contains the
			if (j <= k) l=i;
			//kth element.
		}
	}
	}else{
		printf("SELECT KTH trying find the kth but k>n or n<1 an array with n=%d k=%d \n",n,k);
		exit(1);
	}
}


void *new_Array(int count,size_t element_size)
{
    return calloc(1, sizeof(element_size));
}

void inline quicksortWithWeights(unsigned long n, int arr[],double weights[]){
	if(n>1){
		unsigned long i,ir=n,j,k,l=1,*istack;
		int jstack=0;
		int a,temp;
		istack=lvector(1,NSTACK);
		for (;;) {
			//Insertion sort when subarray small enough.
			if (ir-l < M) {
				for (j=l+1;j<=ir;j++) {
					a=arr[j];
					for (i=j-1;i>=l;i--) {
						if (arr[i] <= a) break;
						arr[i+1]=arr[i];
					}
					arr[i+1]=a;
				}
				if (jstack == 0) break;
				ir=istack[jstack--];
				//Pop stack and begin a new round of parti-
				l=istack[jstack--];
				//tioning.
			} else {
				k=(l+ir) >> 1;
				//Choose median of left, center, and right el-
				SWAPINT(arr[k],arr[l+1]);
				SWAPDOUBLE(weights[k],weights[l+1]);
				//ements as partitioning element a. Also
				if (arr[l] > arr[ir]) {
					//rearrange so that a[l] ≤ a[l+1] ≤ a[ir].
					SWAPINT(arr[l],arr[ir]);
					SWAPDOUBLE(weights[l],weights[ir]);
				}
				if (arr[l+1] > arr[ir]) {
					SWAPINT(arr[l+1],arr[ir]);
					SWAPDOUBLE(weights[l+1],weights[ir]);
				}
				if (arr[l] > arr[l+1]) {
					SWAPINT(arr[l],arr[l+1]);
					SWAPDOUBLE(weights[l],weights[l+1]);
				}
				i=l+1;
				//Initialize pointers for partitioning.
				j=ir;
				a=arr[l+1];
				//Partitioning element.
				for (;;) {
					//Beginning of innermost loop.
					do i++; while (arr[i] < a);
					//Scan up to find element > a.
					do j--; while (arr[j] > a);
					//Scan down to find element < a.
					if (j < i) break;
					//Pointers crossed. Partitioning complete.
					SWAPINT(arr[i],arr[j]);
					SWAPDOUBLE(weights[i],weights[j]);
					//Exchange elements.
				}
				//End of innermost loop.
				arr[l+1]=arr[j];
				//Insert partitioning element.
				arr[j]=a;
				jstack += 2;
				//Push pointers to larger subarray on stack, process smaller subarray immediately.
				if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
				if (ir-i+1 >= j-l) {
					istack[jstack]=ir;
					istack[jstack-1]=i;
					ir=j-1;
				} else {
					istack[jstack]=j-1;
					istack[jstack-1]=l;
					l=i;
				}
			}
		}
		free_lvector(istack,1,NSTACK);
	}else{
		printf(" SORT WITH WEIGHTS trying to short  an array with n=%d \n",n);
		exit -1;
	}
}

int generateRandomPositiveNegPotiveValue(int max,int min){
	if ((max - 1) == RAND_MAX) {
	    return rand();
	  } else {
	    long end = RAND_MAX / max; // truncate skew
	    assert (end > 0L);
	    end *= max;
	    int r;
	    while ((r = rand()) >= end);
	    return r % max;
	  }
}
  int* generateRandomArray(int max,int min,int arraySize){
	int i=0;
	int * arr;
	arr=(int *) malloc (sizeof (int) * arraySize);//new_Array
	for ( i = 0; i < arraySize; i++ )
	{
		 int temp=generateRandomPositiveNegPotiveValue(max,min);
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
		printf("The array is not partitioned well!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		MPI_Scatter(arr, arrayNumberOfPartitionElements, MPI_INT, receiveBuffer,
						arrayNumberOfPartitionElements, MPI_INT, 0, MPI_COMM_WORLD);
	}
}

int findWeightedMedian(int *arr,double *weights,int start,int end){
	if( sizeof(arr) != sizeof(weights) ){
		fprintf(stderr, "Error the arrays must be of equal size \n");
	}else{
		quicksortWithWeights(end,arr-1,weights-1);
		double sum = 0;
		int i;
		for(i = 0; i < end+1; ++i)
		{
			sum += weights[i];
			if(sum >= (double)1/2 )
				break;
		}
		int median = arr[i-1];
		return median;
	}
}

void createMiNiVector(int localNumberOfElements, int kthIndex,
		int* miNiSendBuffer, int* localElementsRcvBuf) {
	if (localNumberOfElements != 0) {
		*(miNiSendBuffer) = localElementsRcvBuf[kthIndex];
		*(miNiSendBuffer + 1) = localNumberOfElements;
	} else {
		*(miNiSendBuffer) = 0;
		*(miNiSendBuffer + 1) = 0;
	}
}

void computeWeights(int p, int numberOfEllementsLeft,
		int* miNiRcvBuffer, double* weights, int* miVector) {
	int ii=0;
	for (ii = 0; ii < p; ii++) {
		int ni = *(miNiRcvBuffer + (ii * 2) + 1);
		int mi = *(miNiRcvBuffer + (ii * 2));
		*(weights + ii) = (double) ni / numberOfEllementsLeft;
		*(miVector + ii) = mi;
		if (ni == 0) {
			*(weights + ii) = 0;
			*(miVector + ii) = 0;
		}
	}
}

LOCALCOUNTS* computeLocalCounts(int localNumberOfElements, int weightedMedian, int* localElementsRcvBuf, LEG* leg)
{
	LOCALCOUNTS* localCounters=calloc(1,sizeof(LOCALCOUNTS));
	localCounters->ei=0;
	localCounters->gi=0;
	localCounters->li=0;
	int ii=0;
	for(ii=0 ; ii< localNumberOfElements ;ii++){
		int localElement=*(localElementsRcvBuf+ii);
		if(localElement < weightedMedian){//less
			*(leg->L+(localCounters->li))=localElement;
			localCounters->li++;
		}else if(localElement > weightedMedian){//greater
			*(leg->G+(localCounters->gi))=localElement;
			localCounters->gi++;
		}else{//equal
			*(leg->E+(localCounters->ei))=localElement;//equa
			localCounters->ei++;
		}
	}
	return localCounters;
}

void createLocalMiNiVector(int localNumberOfElements, int my_rank,
		int* localElementsRcvBuf, int* miNiSendBuffer) {
	if (localNumberOfElements > 1  ) {
		int k= localNumberOfElements / 2;
		int kth = selectkthelem(k,localNumberOfElements,localElementsRcvBuf-1);
		*(miNiSendBuffer) = kth;//after sorting
		*(miNiSendBuffer + 1) = localNumberOfElements;
	} else if (localNumberOfElements==1){
		*(miNiSendBuffer) = *localElementsRcvBuf;
		*(miNiSendBuffer + 1) = 1;
	}
}

int fillBcastBuffer( int p,LOCALCOUNTS *receiveGlobalCounts, int* bcastBuff) {
	int tempCounters[3]={0};
	int ii=0;
	for (ii = 0; ii < p; ii++) {
		tempCounters[0] += ((receiveGlobalCounts+ii)->li);
	}
	*(bcastBuff) = tempCounters[0];
	for (ii = 0; ii < p; ii++) {
		tempCounters[1] += ((receiveGlobalCounts+ii)->gi);
	}
	*(bcastBuff + 1) = tempCounters[1];
	for (ii = 0; ii < p; ii++) {
		tempCounters[2] +=((receiveGlobalCounts+ii)->ei);
	}
	*(bcastBuff + 2) = tempCounters[2];
	return ii;
}
void keepGreaterElements(int localNumberOfElements, int* localElementsRcvBuf,
		LEG* leg) {
	int i = 0;
	for (i = 0; i < localNumberOfElements; i++) {
		*(localElementsRcvBuf + i) = *(leg->G + i);
	}
}
void keepLessElements(int localNumberOfElements, int* localElementsRcvBuf,LEG* leg) {
	int i = 0;
	for (i = 0; i < localNumberOfElements; i++) {
		*(localElementsRcvBuf + i) = *(leg->L + i);
	}
}
  /*
   * input arguments argv[1]=size of random array argv[2]=the rank of the ellement we are looking for
   */
int main(int argc, char* argv[]){
	int my_rank; /* rank of process */
	int p;       /* number of processes */
	int sizeOfMainArray,numberOfEllementsLeft,constantParametrized=100;
	int *localElementsRcvBuf = NULL;
	int *arr = NULL;
	int *miNiRcvBuffer = NULL;
	int *bcastBuff = NULL;
	int *miVector = NULL;
	int *gatherFromAllProcessesArray = NULL;
	int *totalRemainingArrFinal=NULL;
	double *weights = NULL;
	int *miNiSendBuffer = NULL;
	int *elementPerProcessCount=NULL;
	int weightedMedian=0;
	int found=0;
	int intialRankNumber=0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//initialize after the call to MPI_INIT to get the correct argv ....
	sizeOfMainArray=atoi(argv[1]);
	int kSmallestEllement=atoi(argv[2]);
	intialRankNumber=kSmallestEllement;
	int arrayNumberOfPartitionElements=sizeOfMainArray/p;
	LOCALCOUNTS *localCounters=calloc(1,sizeof(LOCALCOUNTS));;
	LOCALCOUNTS *receiveGlobalCounts=calloc(p,sizeof(LOCALCOUNTS));
	LEG *leg=New_LEG(arrayNumberOfPartitionElements);
	NEWARRAY(localElementsRcvBuf, arrayNumberOfPartitionElements);
	NEWARRAY(miVector, p);
	NEWARRAY(bcastBuff, 3);
	NEWARRAY(miNiSendBuffer, 2);
	numberOfEllementsLeft=sizeOfMainArray;//INITIAL VALUE OF REMAINIGN NUMBER OF ELLEMENTS

	/* create a type for struct localCounters */
	const int nitems=3;
	int          blocklengths[3] = {1,1,1};
	MPI_Datatype types[3] = {MPI_INT, MPI_INT,MPI_INT};
	MPI_Datatype mpi_local_counts_type;
	MPI_Aint     offsets[3];
	MPI_Aint 	 addr[3];
	MPI_Get_address(localCounters, &addr[0]);
	MPI_Get_address(&localCounters->li, &addr[1]);
	MPI_Get_address(&localCounters->ei, &addr[2]);
	MPI_Get_address(&localCounters->gi, &addr[3]);
	offsets[0] = addr[1] - addr[0];
	offsets[1] = addr[2] - addr[0];
	offsets[1] = addr[3] - addr[0];
	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_local_counts_type);
	MPI_Type_commit(&mpi_local_counts_type);
	/* create a type for struct localCounters */

	if ( my_rank == 0) {
		NEWARRAY(miNiRcvBuffer,p*2);
		NEWARRAY(arr,sizeOfMainArray);
		arr=generateRandomArray(100000,INT_MIN,sizeOfMainArray);
		NEWARRAY(weights,p);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	scatterTheArray(sizeOfMainArray, p, arrayNumberOfPartitionElements, arr,
			localElementsRcvBuf);
	//print_array(localElementsRcvBuf,arrayNumberOfPartitionElements,my_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	int localNumberOfElements=(arrayNumberOfPartitionElements);
	while( numberOfEllementsLeft >= (sizeOfMainArray/(constantParametrized*p) )){
		createLocalMiNiVector(localNumberOfElements, my_rank,
						localElementsRcvBuf, miNiSendBuffer);
		MPI_Gather(miNiSendBuffer,2,MPI_INT,miNiRcvBuffer,2,MPI_INT,0,MPI_COMM_WORLD);
		if(my_rank==0){
			computeWeights(p, numberOfEllementsLeft,
					miNiRcvBuffer, weights, miVector);
			weightedMedian=findWeightedMedian(miVector,weights,0,p-1);//to check the end size!!!!!!!!!!!!!
		}
		MPI_Barrier(MPI_COMM_WORLD);
		*bcastBuff=weightedMedian;
		MPI_Bcast(bcastBuff,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		weightedMedian=*bcastBuff;
		localCounters=computeLocalCounts(localNumberOfElements, weightedMedian, localElementsRcvBuf,leg);//modifies leg also
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(localCounters,1,mpi_local_counts_type,receiveGlobalCounts,1,mpi_local_counts_type,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		fillBcastBuffer(p,receiveGlobalCounts, bcastBuff);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(bcastBuff,3,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		short caseOne=*bcastBuff < kSmallestEllement && kSmallestEllement <= (*bcastBuff + (*(bcastBuff+2)));
		short caseTwo=*bcastBuff >= kSmallestEllement;
		if(caseOne){
			printf("if %d  \n",my_rank);
			found=1;
			break;
		}else if(caseTwo){
			printf("else f %d\n",my_rank);
			numberOfEllementsLeft=*bcastBuff;
			localNumberOfElements=localCounters->li;
			keepLessElements(localNumberOfElements,
					localElementsRcvBuf,leg);
		}else{
			numberOfEllementsLeft=*(bcastBuff+1);
			localNumberOfElements=localCounters->gi;
			kSmallestEllement=kSmallestEllement-(*bcastBuff+(*(bcastBuff+2)));
			keepGreaterElements(localNumberOfElements,localElementsRcvBuf,leg);
		}
		localCounters->ei=0;
		localCounters->gi=0;
		localCounters->li=0;
		MPI_Barrier(MPI_COMM_WORLD);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	if(!found){
		if(my_rank==0){
			NEWARRAY(gatherFromAllProcessesArray,numberOfEllementsLeft*p);
			NEWARRAY(totalRemainingArrFinal,numberOfEllementsLeft*p);

			NEWARRAY(elementPerProcessCount,p);
		}
		MPI_Gather(&localNumberOfElements,1,MPI_INT,elementPerProcessCount,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(localElementsRcvBuf,numberOfEllementsLeft,MPI_INT,gatherFromAllProcessesArray,numberOfEllementsLeft,MPI_INT,0,MPI_COMM_WORLD);	//here is the mistake
		if(my_rank==0){
			/*here we popullate the array from the different sizes of the received elements from each process into a main array called totalRemainingArrFinal */
			int i=0;
			int jj=0;
			int kk=0;
			while(jj<numberOfEllementsLeft && kk < p ){
				while(i < *(elementPerProcessCount+kk) ){
					*(totalRemainingArrFinal+jj)=*( (gatherFromAllProcessesArray+(kk*numberOfEllementsLeft)+i) );
					jj++;
					i++;
				}
				i=0;
				kk++;
			}
			/*finished populating the final array contains the different sized arrays from all processes
			 *
			 */
			int finalKthElement=selectkthelem(kSmallestEllement,numberOfEllementsLeft,totalRemainingArrFinal-1);
			quicksort(sizeOfMainArray,arr-1);
			quicksort(numberOfEllementsLeft,gatherFromAllProcessesArray-1);
			printf("The number we are looking for is x:=%d and the result from computation is res:=%d \n",arr[intialRankNumber-1],totalRemainingArrFinal[kSmallestEllement-1]);
			printf("the k-th element is %d ",finalKthElement);
		}
	}else{
		if(my_rank==0){
			quicksort(sizeOfMainArray,arr-1);
			printf("The number we are looking for is x:=%d and the result from computation is res:=%d \n",arr[intialRankNumber-1],weightedMedian);
			printf("the k-th element is %d  after tfound the weighted \n",weightedMedian);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
