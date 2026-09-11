#ifndef LPLATTICE_H
#define LPLATTICE_H

/***************** LPLattice ******************
 *  Library to create and manage lattice systems
 *  for simulations (e.g. clusters, percolation,
 *  statistical physics models).
 *
 *  Supports 3D-dimensional lattices with different
 *  boundary conditions and neighbor structures.
 *********************************************/

#include <stdbool.h>
#include <stddef.h>

typedef struct Neighbours {
    size_t size;
    size_t neighboursCount;
    int** matrix;
} Neighbours;

/***************** Dimension *****************
 * Represents the geometry of an N-dimensional lattice.
 *
 * count        -> number of dimensions (e.g. 2D, 3D)
 * size         -> total number of sites in the lattice
 * proportions  -> size per axis (heap-allocated array)
 *
 * Example:
 *   2D lattice 10x20:
 *     count = 2
 *     proportions = {10, 20}
 *     size = 200
 *********************************************/
typedef struct Dimension {
    size_t count;
    size_t size;
    size_t* proportions;
} Dimension;

/***************** BoundaryCondition *****************
 * Defines how lattice edges behave:
 *
 * BC_PERIODIC        -> wrap-around boundaries
 * BC_OPEN_VERTICAL   -> open in vertical direction
 * BC_OPEN_HORIZONTAL -> open in horizontal direction
 * BC_OPEN            -> fully open boundaries
 *****************************************************/
typedef enum {
    BC_PERIODIC = 0,
    BC_OPEN_VERTICAL = 1,
    BC_OPEN_HORIZONTAL = 2,
    BC_OPEN = 3,
} BoundaryCondition;

/***************** Forward declarations *****************
 * Lattice is an opaque type: its internal structure
 * is hidden from users of the API.
 *******************************************************/
typedef struct Lattice {
    Dimension dimension;
    BoundaryCondition boundryCondition;
    Neighbours* neighbours;
} Lattice;

/***************** Dimension API *****************
 * Creates and destroys a Dimension structure.
 *
 * IMPORTANT:
 * - createDimension allocates memory for proportions[]
 * - destroyDimension MUST be called to free it
 *************************************************/
void createDimension(Dimension* dimension, size_t axis_count, const size_t* proportions);
void destroyDimension(Dimension* dimension);

/***************** Lattice API *****************
 * Creates and destroys a lattice system.
 *
 * Lattice owns:
 * - its Dimension (copied by value, but contains heap memory)
 * - its Neighbours structure (internal)
 *************************************************/
void createLattice(Lattice* lattice, const Dimension* dimension,
                   BoundaryCondition boundryCondition);

void destroyLattice(Lattice* lattice);

/***************** Interaction rules *****************
 * Checks whether two sites are valid neighbors
 * under the current boundary condition rules.
 *****************************************************/
bool isValidInteraction(const Lattice* lattice, int s0, int s1);

/***************** Cluster labeling *****************
 * Implements Hoshen–Kopelman algorithm variant.
 *
 * label[i] -> cluster representative of site i
 *****************************************************/
void labelCluster(int* label, const int* s, const Lattice* lattice);

/***************** Debug / output *****************
 * Writes lattice state to file or stdout (implementation-defined).
 *****************************************************/
void plotMatrix(const int* matrix, const Lattice* lattice);

/***************** Neighbour query *****************
 * Returns index of a neighbor of a site.
 *
 * neighbourMatrixIndex mapping:
 *   0 -> RIGHT
 *   1 -> LEFT
 *   2 -> DOWN
 *   3 -> UP
 *****************************************************/
int getNeighbour(const Lattice*, int mainIndex, int neighbourMatrixIndex);

#endif
