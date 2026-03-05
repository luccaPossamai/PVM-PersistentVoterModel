#ifndef LPUTIL_H
#define LPUTIL_H

#include <stddef.h>

unsigned int setupRandom(unsigned int seed);

void lputil_die(const char* errorMessage);

float randomFloat(void);

int randomInt(int endExclusive);

int randomIntOf(int startInclusive, int endExclusive);

void* safeMAlloc(size_t memSize);

int sumIntArray(const int* array, size_t arraySize);

float sumFloatArray(const float* array, size_t arraySize);

void geomProgression(float* array, float startInclusive, float endInclusive, size_t arraySize);

void linearProgression(float* array, float startInclusive, float endInclusive, size_t arraySize);

#endif
