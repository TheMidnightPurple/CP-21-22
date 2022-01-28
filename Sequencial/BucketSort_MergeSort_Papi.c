#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////////////////////
#include "papi.h"

// PAPI events to monitor
#define NUM_EVENTS 4
int Events[NUM_EVENTS] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_L1_DCM, PAPI_L2_DCM};
// PAPI counters' values
long long values[NUM_EVENTS], min_values[NUM_EVENTS];
int retval, EventSet=PAPI_NULL;

// number of times the function is executed and measured
#define NUM_RUNS 5
////////////////////////////////////////////////////////////////////////////////////////

#define MAX_RANGE 100000
#define MIN_RANGE 1

typedef struct {
    int actual_column_index;
    int number_elements;
    int* inside_bucket;
} Bucket;

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    int* L = malloc(n1*sizeof(int));
    int* R = malloc(n2*sizeof(int));

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
            if (L[i] <= R[j]) {
                    arr[k] = L[i];
                    i++;
            }
            else {
                    arr[k] = R[j];
                    j++;
            }
            k++;
    }

    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
            arr[k] = L[i];
            i++;
            k++;
    }

    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
            arr[k] = R[j];
            j++;
            k++;
    }

    free(L);
    free(R);
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
  
        merge(arr, l, m, r);
    }
}

void bucketSort(Bucket* bucket, int* vector, int n_buckets, int vector_size, int max_value){
    //Ordenar cada Bucket e escrever resultado num ficheiro(para ter paralelismo no radixSort, tinha de separar o sort da escrita no ficheiro)
    int count = 0;

    for (int i=0; i<n_buckets; i++){
        if (bucket[i].number_elements != 0){
            mergeSort(bucket[i].inside_bucket, 0, bucket[i].number_elements - 1);
        }
    }

    for (int i=0; i<n_buckets; i++){
        for(int j=0; j < bucket[i].number_elements; j++){
            vector[count++] = bucket[i].inside_bucket[j];
        }
    }
}

void bucketInit(Bucket* bucket, int* vector, int n_buckets, int vector_size, int max_value){

    //Inicilizar coluna atual e numero de elementos em todos os buckets
    for(int i=0; i<n_buckets; i++){
        bucket[i].actual_column_index = 0;
        bucket[i].number_elements = 0;
    }
}

int generateValues(int* vector, int numElements){

    int max_value = 0;

    for (int i=0; i<numElements; i++){
        int number = rand() % MAX_RANGE + MIN_RANGE;
        vector[i] = number;

        if (number >= max_value){
            max_value = number;
        }
    }

    return max_value;
}

void bucketPutValue(Bucket* bucket, int* vector, int n_buckets, int vector_size, int max_value){
    //Colocar cada valor no seu Bucket(pode ter paralelismo)
    for (int i=0; i<vector_size; i++){
        int index = floor(n_buckets * vector[i] / max_value);

        if (index >= (n_buckets)){
            index--;//para não haver casos onde passa do limite
        }

        bucket[index].number_elements++;
    }

    for (int i=0; i<n_buckets; i++){

        int n_elements = bucket[i].number_elements;
        bucket[i].inside_bucket = malloc(n_elements * sizeof(int));
    }

    for (int i=0; i<vector_size; i++){

        int index = floor(n_buckets * vector[i] / max_value);

        if (index >= (n_buckets)){
            index--;//para não haver casos onde passa do limite
        }

        int column = bucket[index].actual_column_index;
        bucket[index].inside_bucket[column] = vector[i];
        bucket[index].actual_column_index++;
    }
}


int mainPAPI(int NUMBER_BUCKETS, int NUMBER_ELEMENTS){

    long long start_usec, end_usec, elapsed_usec, min_usec=0L;
    int i, run;
    int num_hwcntrs = 0;

    fprintf (stdout, "Array have a total of %d elements!\n", NUMBER_ELEMENTS);

    fprintf (stdout, "\nSetting up PAPI...");
    // Initialize PAPI
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fprintf(stderr,"PAPI library init error!\n");
        return 0;
    }

    /* create event set */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr,"PAPI create event set error\n");
        return 0;
    }


    /* Get the number of hardware counters available */
    if ((num_hwcntrs = PAPI_num_hwctrs()) <= PAPI_OK)  {
        fprintf (stderr, "PAPI error getting number of available hardware counters!\n");
        return 0;
    }
    fprintf(stdout, "done!\nThis system has %d available counters.\n\n", num_hwcntrs);

    // We will be using at most NUM_EVENTS counters
    if (num_hwcntrs >= NUM_EVENTS) {
        num_hwcntrs = NUM_EVENTS;
    } else {
        fprintf (stderr, "Error: there aren't enough counters to monitor %d events!\n", NUM_EVENTS);
        return 0;
    }

    if (PAPI_add_events(EventSet,Events,NUM_EVENTS) != PAPI_OK)  {
        fprintf(stderr,"PAPI library add events error!\n");
        return 0;
    }

    // init array
    fprintf (stdout, "Initializing Array with random values...");
    int max_value = 0;
    int* vector = malloc(NUMBER_ELEMENTS * sizeof(int));
    max_value = generateValues(vector, NUMBER_ELEMENTS);
    fprintf (stdout, "done!\n");

    // ini buckets
    fprintf (stdout, "Create Bucket...");
    Bucket* bucket = malloc(NUMBER_BUCKETS * sizeof(Bucket));
    fprintf (stdout, "done!\n");

    // warmup caches
    fprintf (stdout, "Warming up caches...");
    bucketInit(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    bucketPutValue(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    bucketSort(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    fprintf (stdout, "done!\n");


    for (run=0 ; run < NUM_RUNS ; run++) {
        fprintf (stdout, "run=%d - Sorting Array...", run);

        // use PAPI timer (usecs) - note that this is wall clock time
        // for process time running in user mode -> PAPI_get_virt_usec()
        // real and virtual clock cycles can also be read using the equivalent
        // PAPI_get[real|virt]_cyc()
        start_usec = PAPI_get_real_usec();

        /* Start counting events */
        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf (stderr, "PAPI error starting counters!\n");
            return 0;
        }

        bucketInit(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
        bucketPutValue(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
        bucketSort(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);

        /* Stop counting events */
        if (PAPI_stop(EventSet,values) != PAPI_OK) {
            fprintf (stderr, "PAPI error stoping counters!\n");
            return 0;
        }

        end_usec = PAPI_get_real_usec();
        fprintf (stdout, "done!\n");

        elapsed_usec = end_usec - start_usec;

        if ((run==0) || (elapsed_usec < min_usec)) {
            min_usec = elapsed_usec;
            for (i=0 ; i< NUM_EVENTS ; i++) min_values[i] = values [i];
        }

    } // end runs
    fprintf (stdout,"\nWall clock time: %lld usecs\n", min_usec);

    // output PAPI counters' values
    for (i=0 ; i< NUM_EVENTS ; i++) {
        char EventCodeStr[PAPI_MAX_STR_LEN];

        if (PAPI_event_code_to_name(Events[i], EventCodeStr) == PAPI_OK) {
            fprintf (stdout, "%s = %lld\n", EventCodeStr, min_values[i]);
        } else {
            fprintf (stdout, "PAPI UNKNOWN EVENT = %lld\n", min_values[i]);
        }
    }

#if NUM_EVENTS >1
  // evaluate CPI and Texec here
  if ((Events[0] == PAPI_TOT_CYC) && (Events[1] == PAPI_TOT_INS)) {
    float CPI = ((float) min_values[0]) / ((float) min_values[1]);
    fprintf (stdout, "CPI = %.2f\n", CPI);
  }
#endif

    /*for (int i=0; i<NUMBER_ELEMENTS; i++){
        fprintf (stdout, "%d\n", vector[i]);
    }*/

    fprintf (stdout,"\nThat's all, folks\n");

    free(vector);

    for (int i = 0; i < NUMBER_BUCKETS; ++i)
    {
        free(bucket[i].inside_bucket);
    }
    free(bucket);

    return 0;
}

int mainSemPAPI(int NUMBER_BUCKETS, int NUMBER_ELEMENTS){
    

    // init array
    fprintf (stdout, "Initializing Array with random values...");
    int max_value = 0;
    int* vector = malloc(NUMBER_ELEMENTS * sizeof(int));
    max_value = generateValues(vector, NUMBER_ELEMENTS);
    fprintf (stdout, "done!\n");

    // ini buckets
    fprintf (stdout, "Create Bucket...");
    Bucket* bucket = malloc(NUMBER_BUCKETS * sizeof(Bucket));
    fprintf (stdout, "done!\n");

    // warmup caches
    fprintf (stdout, "Playing...");
    
    bucketInit(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    bucketPutValue(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    bucketSort(bucket, vector, NUMBER_BUCKETS, NUMBER_ELEMENTS, max_value);
    
    fprintf (stdout, "done!\n");

    /*for (int i=0; i<NUMBER_ELEMENTS; i++){
        fprintf (stdout, "%d\n", vector[i]);
    }*/

    free(vector);

    for (int i = 0; i < NUMBER_BUCKETS; ++i)
    {
        free(bucket[i].inside_bucket);
    }
    free(bucket);

    return 0;
}


int verify_command_line (int argc, char *argv[], int *total_elements, int *buckets, int *tipo) {
    int val, val2, val3;

    if (argc!=4) {
            fprintf (stdout, "Exactly 3 arguments are required!");
            return 0;
    }


    val = atoi (argv[1]);
    if (val <= 0) {
            fprintf (stdout, "The vector size is the size of vector and must be a positive integer!");
            return 0;
    }
    else {
            *total_elements = val;
    }


    val2 = atoi (argv[2]);
    if (val2 <= 0) {
            fprintf (stdout, "The number of buckets are the number of buckets used and must be a positive integer!");
            return 0;
    }
    else {
            *buckets = val2;
    }


    val3 = atoi (argv[3]);
    if (val3 != 0 && val3 != 1) {
            fprintf (stdout, "The number of type must be 0 or 1! (0 = use PAPI | 1 = without PAPI)");
            return 0;
    }
    else {
            *tipo = val3;
    }

    return 1;
}


int main(int argc, char *argv[]){
    
    int NUMBER_BUCKETS;
    int NUMBER_ELEMENTS;
    int TIPO;

    if (!verify_command_line (argc, argv, &NUMBER_ELEMENTS, &NUMBER_BUCKETS, &TIPO)) {
            return 1;
    }


    if(TIPO == 0){
        mainPAPI(NUMBER_BUCKETS, NUMBER_ELEMENTS);
    }else{
        if(TIPO == 1){
            mainSemPAPI(NUMBER_BUCKETS, NUMBER_ELEMENTS);
        }
    }
}