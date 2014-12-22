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

#define SWAPINT(r,s)  do{int t=r; r=s; s=t; } while(0)
#define SWAPDOUBLE(r,s)  do{double t=r; r=s; s=t; } while(0)
#define NEWARRAY(x, n) do { (x) = calloc((n), sizeof *(x)); } while (0)
#include <time.h>


void *new_Array(int count,size_t element_size)
{
    return calloc(1, sizeof(element_size));
error:
    return NULL;
}

void quicksort(int *array, int start, int end){
    if(start < end){
        int l=start+1, r=end, p = array[start];
        while(l<r){
            if(array[l] <= p)
                l++;
            else if(array[r] >= p)
                r--;
            else
            	SWAPINT(array[l],array[r]);
        }
        if(array[l] < p){
        	SWAPINT(array[l],array[start]);
            l--;
        }
        else{
            l--;
            SWAPINT(array[l],array[start]);
        }
        quicksort(array, start, l);
        quicksort(array, r, end);
    }
}
void quicksortWithWeights(int *array,double *weights, int start, int end){
    if(start < end){
        int l=start+1, r=end, p = array[start];
        while(l<r){
            if(array[l] <= p)
                l++;
            else if(array[r] >= p)
                r--;
            else
            	SWAPINT(array[l],array[r]);
            	SWAPDOUBLE(weights[l],weights[r]);
        }
        if(array[l] < p){
        	SWAPINT(array[l],array[start]);
        	SWAPDOUBLE(weights[l],weights[start]);
            l--;
        }
        else{
            l--;
            SWAPINT(array[l],array[start]);
            SWAPDOUBLE(weights[l],weights[start]);
        }
        quicksortWithWeights(array,weights, start, l);
        quicksortWithWeights(array,weights, r, end);
    }
}

// selects the median of medians in an array
int selectkth(int *a, int s, int e, int k,int my_rank){
    // if the partition length is less than or equal to 5
    // we can sort and find the kth element of it
    // this way we can find the median of n/5 partitions
	//printf("endering select s=%d k=%d e=%d my_rank=%d \n\n",s,k,e,my_rank);
    if(e-s+1 <= 5){
    	int i=0;
    	for(i=0; i<e+1 ; i++ ){
    	    		//printf("selectkth before sort i=%d value=%d my_rank=%d \n",i,a[i],my_rank);
    	    	}
    	//printf("\n\n");
    	quicksort(a,s,e-1);
    	 i=0;
    	for( i=0;i<e+1;i++ ){
    		//printf("selectkth after sort i=%d value=%d my_rank=%d \n",i,a[i],my_rank);
    	}
    	//printf("returning from here s=%d k=%d my_rank=%d \n\n",s,k,my_rank);
    	if(k+s==0){
    		return 0;
    	}else{
    		return s+k-1;
    	}
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
        int median = selectkth(a, 5*i, 5*i+4, 3,my_rank);
        SWAPINT(a[median], a[i]);
    }

    // now we have array
    // a[0] = median of 1st 5 sized partition
    // a[1] = median of 2nd 5 sized partition
    // and so on till n/5
    // to find out the median of these n/5 medians
    // we need to select the n/10th element of this set (i.e. middle of it)
    return selectkth(a, 0, (e+1)/5, (e+1)/10,my_rank);
}

int generateRandomPositiveNegPotiveValue(int max,int min){
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
	arr=(int *) malloc (sizeof (int) * arraySize);//new_Array
	for ( i = 0; i < arraySize; i++ )
	{
		 int temp=generateRandomPositiveNegPotiveValue(max,min);
		// printf("random is %d:=%d \n",i,temp);
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

void printEllementsOfProsseses(int my_rank, int arrayNumberOfPartitionElements,
		int* receiveBuffer) {
	//struct timespec tim, tim2;
	//tim.tv_sec = 1;
		//tim.tv_nsec = 5000;

	if (my_rank != 0) {
		//int ii = 0;
		//for (ii = 0; ii < arrayNumberOfPartitionElements; ii++) {
			//printf("element %d=%d process %d ", ii, (receiveBuffer[ii]),
			//		my_rank);
			//printf("\n");
		//}
		//printf("\n");
	} else {
		//int ii = 0;
		//printf("Process %d received data \n", my_rank);
		//for (ii = 0; ii < arrayNumberOfPartitionElements; ii++) {
			//printf("element %d=%d process %d ", ii, (receiveBuffer[ii]),
			//		my_rank);
			//printf("\n");
		//}
		//printf("\n");
	}

}

int findWeightedMedian(int *arr,double *weights,int start,int end){
	if( sizeof(arr) != sizeof(weights) ){
		fprintf(stderr, "Error the arrays must be of equal size \n");
	}else{
		int ii = 0;
		/*
		for (ii = 0; ii < end+1; ii++) {
					printf("findWeightedMedian element %d=%f process %d ", ii, (weights[ii]),
							0);
					printf("\n");
		}
		printf("\n");
		*/
		quicksortWithWeights(arr,weights,start,end);
		/*
		for (ii = 0; ii < end+1; ii++) {
							printf("findWeightedMedia end=%d AFTER QUICK SORT element %d=%f process %d ",end, ii, (weights[ii]),
									0);
							printf("\n");
				}

		printf("\n");

		for (ii = 0; ii < end+1; ii++) {
									printf("findWeightedMedia end=%d AFTER QUICK SORT element %d=%f process %d ",end, ii, (arr[ii]),
											0);
									printf("\n");
						}

				printf("\n");
				*/


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

void print_array(int *array, int length,int rank)
{
	//int i;
    //for ( i = 0; i < length; i++) { printf("element %d=%d process %d \n", i ,array[i] ,rank); }
    //printf("\n");
}
void print_array_global(int *array, int length,int rank,int true){
	int i;
	if(true){
	    //for ( i = 0; i < length; i++) { printf("GLOBAL element %d=%d process %d \n", i ,array[i] ,rank); }
	    //printf("\n");
	}
	else{
		//for ( i = 0; i < length; i++) { printf("THE BCAST BUFFER  element %d=%d process %d \n", i ,array[i] ,rank); }
		//printf("\n");
	}
}

  /*
   * input arguments argv[1]=size of random array argv[2]=the rank of the ellement we are looking for
   */
int main(int argc, char* argv[]){
	 clock_t start, end,start1,end1,start2,end2,start3,end3,start4,end4;
	 double cpu_time_used,cpu_time_used1;
	 start = clock();


	int my_rank; /* rank of process */
	int p;       /* number of processes */
	int sizeOfMainArray,numberOfEllementsLeft,constantParametrized=100;
	int *localElementsRcvBuf, *arr, *miNiRcvBuffer, *bcastBuff, *L,*E,*G,*Lreceive,*Ereceive,*Greceive,*miVector,*totalRemainingArr;
	double *weights;
	int weightedMedian=0;
	int found=0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//initialiate after the call to MPI_INIT to get the correct argv ....
	struct timespec tim, tim2;
	tim.tv_sec = 1;
	tim.tv_nsec = 5000;

	sizeOfMainArray=atoi(argv[1]);
	int arrayNumberOfPartitionElements=sizeOfMainArray/p;
	//l_ele_rcv_buf=(int *)malloc(arrayNumberOfPartitionElements*sizeof(int));/// all allocate a place holder for scattered data
	NEWARRAY(localElementsRcvBuf,arrayNumberOfPartitionElements);
	NEWARRAY(L,arrayNumberOfPartitionElements);
	NEWARRAY(E,arrayNumberOfPartitionElements);
	NEWARRAY(G,arrayNumberOfPartitionElements);
	NEWARRAY(Lreceive,1);
	NEWARRAY(Ereceive,1);
	NEWARRAY(Greceive,1);
	NEWARRAY(miVector,p);
	//bcast_buff=(int *)malloc(sizeof(int));
	NEWARRAY(bcastBuff,3);
	numberOfEllementsLeft=sizeOfMainArray;//INITIAL VALUE OF REMAINIGN NUMBER OF ELLEMENTS
	if ( my_rank == 0) {
	       //mi_ni_rcv_buff=(int*)malloc(p*2*sizeof(int));
	       NEWARRAY(miNiRcvBuffer,p*2);
	}

	int kSmallestEllement=atoi(argv[2]);
	if (my_rank !=0){//slaves ---

	}
	else{//master
		NEWARRAY(arr,sizeOfMainArray);
		arr=generateRandomArray(100000,INT_MIN,sizeOfMainArray);
		//printf("Created random array \n");
		NEWARRAY(weights,p);
		int i;
		//for ( i = 0; i < sizeOfMainArray; i++) { printf("MAIN ARRAY %d=%d process %d \n", i ,arr[i] ,0); }
	}
	int *miNiSendBuffer;
	NEWARRAY(miNiSendBuffer,2);///every loop allocate new one !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	scatterTheArray(sizeOfMainArray, p, arrayNumberOfPartitionElements, arr,
			localElementsRcvBuf);
	//printf("scattered the array \n\n");

	MPI_Barrier(MPI_COMM_WORLD);
	//printEllementsOfProsseses(my_rank, arrayNumberOfPartitionElements,
			//localElementsRcvBuf);//tracing
	int localNumberOfElements=(arrayNumberOfPartitionElements);

	//printEllementsOfProsseses(my_rank,arrayNumberOfPartitionElements,localElementsRcvBuf);


	while( numberOfEllementsLeft >= (sizeOfMainArray/(constantParametrized*p) )){

		//printf("TestSTARTLOOP");
		//printf("START WHILE numberOfEllementsLeft=%d my_rank%d \n",numberOfEllementsLeft,my_rank);
		//printf("START WHILE localNumberOfElements=%d my_rank%d \n",localNumberOfElements,my_rank);
		start1 = clock();

		int tempcounter=0;
		MPI_Barrier(MPI_COMM_WORLD);



		if(localNumberOfElements!=0){
			/*
			int i=0;
			for ( i = 0; i < localNumberOfElements; i++) { printf("before select kth EIMAI PAPARAS element %d=%d process %d \n", i ,localElementsRcvBuf[i] ,my_rank); }
						    printf("\n");
			*/
			if((localNumberOfElements==1)||(localNumberOfElements/2)==2){
				//printf("INSIDE THE IFIFIFIIFIFIFIFI =%d my_rank=%d \n",localNumberOfElements,my_rank);
			}
			//printf("THE LOCAL NUMBE OF ELEMENTS=%d my_rank=%d \n",localNumberOfElements,my_rank);
			int kthIndex=selectkth(localElementsRcvBuf,0,localNumberOfElements+1,localNumberOfElements/2,my_rank);
			//printf("SELECT KTH RETURED KTH INDEX =%d my_rank=%d \n",kthIndex,my_rank);
			if(localElementsRcvBuf[kthIndex]==0){
				int i=0;
				for(i=0;i<localNumberOfElements;i++)
				printf("the element i=%d value=%d myrank=%d",i,localElementsRcvBuf[i],my_rank);
			}
			printf("\n");
			*(miNiSendBuffer)=localElementsRcvBuf[kthIndex];
			//printf("LOCAL NUMBER OF ELEMENTS =%d my_rank=%d \n",localNumberOfElements,my_rank);
			//printf("KTH INNDEX =%d my_rank=%d \n",kthIndex,my_rank);
			//printf("SELECT KTH RETURED =%d my_rank=%d \n",*miNiSendBuffer,my_rank);
			*(miNiSendBuffer+1)=localNumberOfElements;

			print_array(localElementsRcvBuf,localNumberOfElements,my_rank);
			/*
			for ( i = 0; i < localNumberOfElements; i++) { printf("AFTER   kths EIMAI PAPARAS element %d=%d process %d \n", i ,localElementsRcvBuf[i] ,my_rank); }
			    printf("\n");
			*/

		}else{

			//printf("in here LOCAL NUMBER OF ELEMENTS =0 my_rak%d \n",my_rank);
			*(miNiSendBuffer)=0;
			*(miNiSendBuffer+1)=0;
			//TO DO if this process has o elements
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(miNiSendBuffer,2,MPI_INT,miNiRcvBuffer,2,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		//printf("AFTER THE GATHER rank %d \n",my_rank);
		if(my_rank==0){
			for(tempcounter =0 ;tempcounter <p ;tempcounter++){

				int ni= *(miNiRcvBuffer+(tempcounter*2)+1);
				int mi=*(miNiRcvBuffer+(tempcounter*2));

				//printf(" ni=%d---> mi=%d \n",ni,mi);

				*(weights+tempcounter)= (double) ni / numberOfEllementsLeft;
				*(miVector+tempcounter)= mi;

				if(ni==0){
					*(weights+tempcounter)= (double) ni / numberOfEllementsLeft;
					*(miVector+tempcounter)= mi;
					//printf("in here LOCAL NUMBER OF ELEMENTS =0 my_rak%d \n",my_rank);////////////////////////////////
				}

			}
			int i=0;
			for(i=0;i<p;i++){
							//printf("the mi vector i=%d=%d \n",i,miVector[i]);
							//printf("the weights vector i=%d=%f \n",i,weights[i]);
			}
			weightedMedian=findWeightedMedian(miVector,weights,0,p-1);//to check the end size!!!!!!!!!!!!!

			for(i=0;i<p;i++){
				//printf("the mi vector i=%d=%d \n",i,*miVector);
			}
			//printf("the weightedMedian1 =%d \n",weightedMedian);
		}




		MPI_Barrier(MPI_COMM_WORLD);


		*bcastBuff=weightedMedian;


		//printf("the bcast buff before =%d from rank %d \n",*bcastBuff,my_rank);


		MPI_Bcast(bcastBuff,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

		weightedMedian=*bcastBuff;

		//printf("the weighted median buff =%d \n",weightedMedian);

		int localCounters[3]={0};//l g e counters
		tempcounter=0;
		for(tempcounter=0 ; tempcounter< localNumberOfElements ;tempcounter++){
			int localElement=*(localElementsRcvBuf+tempcounter);
			if(localElement < weightedMedian){//less
				*(L+localCounters[0])=localElement;
				localCounters[0]++;
			}else if(localElement > weightedMedian){//greater
				*(G+localCounters[1])=localElement;
				localCounters[1]++;
			}else{//equal
				*(E+localCounters[2])=localElement;//equa
				localCounters[2]++;
			}
		}


		//print_array(localCounters,3,my_rank);////first PRINT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(&localCounters[0],1,MPI_INT,Lreceive,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(&localCounters[1],1,MPI_INT,Greceive,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(&localCounters[2],1,MPI_INT,Ereceive,1,MPI_INT,0,MPI_COMM_WORLD);
		/*
		if(my_rank==0){
			for(tempcounter=0;tempcounter<p;tempcounter++){
				printf("leceive is =%d FROM RANK %d \n",*(Lreceive+tempcounter),my_rank);////first PRINT!!!!!!
				printf("Greceive IS =%d FROM RANK %d \n",*(Greceive+tempcounter),my_rank);////first PRINT!!!!!!
				printf("Ereceive  N IS =%d FROM RANK %d \n",*(Ereceive+tempcounter),my_rank);////first PRINT!!!!!!
			}
		}
		*/


		tempcounter=0;
		int tempCounters[3]={0};
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank==0){
			for( tempcounter=0 ; tempcounter<p ;tempcounter++){
				tempCounters[0]+=*(Lreceive+tempcounter);
			}
			*(bcastBuff)=tempCounters[0];
			for( tempcounter=0 ; tempcounter<p ;tempcounter++){
				tempCounters[1]+=*(Greceive+tempcounter);
			}
			*(bcastBuff+1)=tempCounters[1];
			for( tempcounter=0 ; tempcounter<p ;tempcounter++){
				tempCounters[2]+=*(Ereceive+tempcounter);
			}
			*(bcastBuff+2)=tempCounters[2];
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//printf("THE REMAINING  N IS =%d FROM RANK %d \n",numberOfEllementsLeft,my_rank);////first PRINT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		MPI_Bcast(bcastBuff,3,MPI_INT,0,MPI_COMM_WORLD);
		//printf("THE REMAINING NUMBERS=%d ",numberOfEllementsLeft);
		print_array_global(bcastBuff,3,my_rank,0);
		int i=0;
		if(*bcastBuff < kSmallestEllement && kSmallestEllement <= (*bcastBuff + (*(bcastBuff+2)))){
			printf("if %d  \n",my_rank);
			found=1;
			break;
		}else if( *bcastBuff >= kSmallestEllement ){
			printf("else f %d\n",my_rank);
			numberOfEllementsLeft=*bcastBuff;
			localNumberOfElements=localCounters[0];
			int i=0;
			for(i=0;i<localNumberOfElements;i++){
				*(localElementsRcvBuf+i)=*(L+i);
				//printf("my elements are i=%d VALUES=%d myrank=%d",i,localElementsRcvBuf[i],my_rank);
			}
		}else{
			//printf("is k>L+E? %d  \n",kSmallestEllement > bcastBuff[0]+bcastBuff[2]);
			printf("else %d  \n",my_rank);
			numberOfEllementsLeft=*(bcastBuff+1);
			localNumberOfElements=localCounters[1];
			kSmallestEllement=kSmallestEllement-(*bcastBuff+*(bcastBuff+2));

			int i=0;
			for(i=0;i<localNumberOfElements;i++){
				*(localElementsRcvBuf+i)=*(G+i);
			}
		}
		end1 = clock();
		memset(localCounters, 0, 3);//reset to zero
		MPI_Barrier(MPI_COMM_WORLD);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("I HAVE EXITED THE LOOP %d",my_rank);

	if(!found){
		if(my_rank==0){
			NEWARRAY(totalRemainingArr,numberOfEllementsLeft*p);
		}

		MPI_Gather(localElementsRcvBuf,numberOfEllementsLeft,MPI_INT,totalRemainingArr,numberOfEllementsLeft,MPI_INT,0,MPI_COMM_WORLD);	//
		if(my_rank==0){
			quicksort(arr,0,sizeOfMainArray-1);
			print_array(arr,sizeOfMainArray,0);
			printf("The 700 smallest element=%d \n",arr[700]);
			printf("The 699 smallest element=%d \n",arr[699]);
			printf("The 701 smallest element=%d \n",arr[701]);
			int finalKthElement=selectkth(totalRemainingArr,0,numberOfEllementsLeft,kSmallestEllement+1,0);
			printf("the k-th element is %d ",totalRemainingArr[finalKthElement]);
		}
	}else{
		if(my_rank==0){
			quicksort(arr,0,sizeOfMainArray-1);
			printf("The 700 smallest element=%d \n",arr[700]);
			printf("The 699 smallest element=%d \n",arr[699]);
			printf("The 701 smallest element=%d \n",arr[701]);
			printf("the k-th element is %d  after tfound the weighted",weightedMedian);
		}
	}
	if(my_rank==0){
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("the time to select %f",cpu_time_used);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	
	return 0;
}
