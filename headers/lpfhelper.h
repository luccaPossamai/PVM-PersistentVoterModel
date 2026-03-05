#ifndef LPFHELPER_H
#define LPFHELPER_H

#include <stdio.h>
#include <stdbool.h>

FILE* safeOpen(const char* filename,
               const char* fileExtension);

FILE* safeSeedOpen(const char* filename,
                   const char* fileExtension,
                   unsigned int* seed,
                   bool forceSeed);

void writeHeader(FILE* file,
                 const char* titleString,
                 const char* authorString);

#endif
