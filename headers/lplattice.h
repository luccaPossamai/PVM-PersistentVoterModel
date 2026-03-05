#ifndef LPLATTICE_H
#define LPLATTICE_H

/*****************LPLattice******************
    Create lattice systems more easily
********************************************/

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    BC_PERIODIC = 0,
    BC_OPEN_VERTICAL = 1,
    BC_OPEN_HORIZONTAL = 2,
    BC_OPEN = 3,
} BoundaryCondition;

/*****************Neighbours*****************
 *     Each agent has 'neighboursCount' neighbours
 *********************************************/
typedef struct Neighbours {
    size_t size;
    size_t neighboursCount;
    int** matrix;
} Neighbours;

/****************SquareLattice***************
    A square lattice can have different dimensionality
    The boundryCondition reflects the neighbours
    some agent has. L is the unidimensional measure
    of the square dimensional lattice.
********************************************/
typedef struct SquareLattice {
    size_t dimension;
    BoundaryCondition boundryCondition;
    size_t L;
    size_t size;
    Neighbours* neighbours;
} SquareLattice;

void initiateSquareLattice(SquareLattice* lattice, size_t dimension, size_t length,
                           BoundaryCondition boundryCondition);
void clearSquareLattice(SquareLattice* lattice);
bool isValidInteraction(const SquareLattice* lattice, int s0, int s1);
void labelCluster(int* labelArray, const int* array, const SquareLattice* lattice);
void labelClusterEff(int* label, const int* s, const SquareLattice* lattice);
void plotMatrix(const int* matrix, const SquareLattice* lattice);

// DONT USE THESE
int getInterfacialLengthWith(int*, SquareLattice*, int, int);
int getInterfacialLength(int*, SquareLattice*, int);
void clusterSizeInfo(int*, int*, int*, int*, int*, SquareLattice*);

#endif
