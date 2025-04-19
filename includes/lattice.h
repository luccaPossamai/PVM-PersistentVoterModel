#include <monte_carlo.h>
#include <math.h>

typedef struct {
    int size;
    int neighboursCount;
    int **matrix;
} Neighbours;

typedef struct {
    int dimension;  
    int L;
    int size;       // elements size
    Neighbours neighbours;     // neighbours matrix
} SquareLattice;



void initiateSquareLattice(SquareLattice*, int, int);
void initiateNeighbours(Neighbours*, int, int);
int validInteraction(SquareLattice*, int, int, int);

void initiateSquareLattice(SquareLattice *l, int dimension, int L){
    int size = pow(L, dimension);
    Neighbours neig;
    initiateNeighbours(&neig, dimension, L);
    l -> dimension = dimension;
    l -> size = size;
    l -> neighbours = neig; 
}
void initiateNeighbours(Neighbours *neig, int dimension, int L){
    neig -> size = pow(L, dimension);
    neig -> neighboursCount = pow(2, dimension);
    neig -> matrix = smalloc(neig -> size * sizeof(int*));
    
    for(int i = 0; i < neig->size; i++){
        neig -> matrix[i] = smalloc(neig->neighboursCount * sizeof(int));
    }
    for(int i = 0; i < neig->size; i++){
        neig -> matrix[i][0] = i % L == (L - 1) ? i - L + 1 : i + 1;
        neig -> matrix[i][1] = i % L == 0       ? i + L - 1 : i - 1;
        if(dimension > 1){
            int N = pow(L, 2);
            neig -> matrix[i][2] = i < L        ? i + N - L : i - L;
            neig -> matrix[i][3] = i >= (N - L) ? i - N + L : i + L;  
        }
    }
}
int validInteraction(SquareLattice* lattice, int s1, int s2, int b){
    if(lattice->dimension == 1){
        int f = abs(s1 - s2) > 1; // means its from periodic sides;
        if(f && b == 3){
            return 0;
        }
    
    } else {
        
    }

    return 1;
}


