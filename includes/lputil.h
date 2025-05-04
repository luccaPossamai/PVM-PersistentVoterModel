#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

unsigned int setupRandom(unsigned int);
float randomFloat(void);
int randomInt(int);
int randomIntOf(int, int);
void *safeMAlloc(unsigned int);
void safeMFree(void*);
int sumIntArray(int*, int);
float sumFloatArray(float*, int);
void geomProgression(float*, float, float, int);
void linearProgression(float*, float, float, int);

unsigned int setupRandom(unsigned int seed){
  if(seed == 0){
    seed = time(NULL);
  }
	srand(seed);
  return seed;
}

float randomFloat(void){
  return (float)rand()/((float)RAND_MAX);
}
int randomInt(int endExclusive){
  return randomIntOf(0, endExclusive);
}
int randomIntOf(int startInclusive, int endExclusive){
  int size = endExclusive - startInclusive;
  return startInclusive + (rand() % size);
}



/************************************************
		Safe Memory Allocation
			-lpcamors(11/24)
	Triggers an error message if malloc(int) fails  
 ************************************************/
void *safeMAlloc(unsigned int size){
  void *ptr = malloc(size);
  if (ptr == NULL) {
      printf("Error: malloc got NULL!\n");
      exit(EXIT_FAILURE);
  }
  return ptr;
}

void safeMFree(void* ptr){
    if(ptr == NULL){
        printf("Error: free got NULL as a parameter!\n");
        return;
    }
    free(ptr);
}


/************************************************
		Sum of Integers of an array
			-lpcamors(11/24)  
 ************************************************/
int sumIntArray(int *arr, int size){
  int sum = 0;
  for(int i = 0; i < size; i ++){
    sum += arr[i];
  }
  return sum;
}
/************************************************
		Sum of Floats of an array
			-lpcamors(11/24)  
 ************************************************/
float sumFloatArray(float *arr, int size){
  float sum = 0;
  for(int i = 0; i < size; i ++){
    sum += arr[i];
  }
  return sum;
}

/************************************************
		Geometrical Progression Generator
			-lpcamors(11/24)
	Maps the "*arr" with geometrically distributed
	elements between "start" and "end" (inclusive)  
 ************************************************/
void geomProgression(float *arr, float start, float end, int size){
  int N = size - 1;
  float a = pow(end / start, 1 / (float) N);// 2^{size + 1} = x
  for(int i = 0; i < size; i++){
    arr[i] = start * pow(a, i);
  }
  
}

/************************************************
		Linear Progression Generator
			-lpcamors(11/24)
	Maps the "*arr" with linear distributed
	elements between "start" and "end" (inclusive)  
 ************************************************/
void linearProgression(float *arr, float start, float end, int size){
  int N = size - 1;
  float a = (end - start) / (float) N;
  for(int i = 0; i < size; i++){
    arr[i] = start + (a * (float) i);
  }
}

