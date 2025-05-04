#include <lputil.h>
#include <math.h>

/*****************LPLattice******************
    Create lattice systems more easily
********************************************/

/*****************Neighbours*****************
    Each agent has 'neighboursCount' neighbours
********************************************/
typedef struct {
    int size;
    int neighboursCount;
    int **matrix;
} Neighbours;



/****************SquareLattice***************
    A square lattice can have different dimensionality
    The boundryCondition reflects the neighbours
    some agent has. L is the unidimensional measure
    of the square dimensional lattice.
********************************************/
typedef struct {
    int dimension;
    int boundryCondition;  
    int L;
    int size;       // elements size
    Neighbours* neighbours;     // neighbours matrix
} SquareLattice;



void initiateSquareLattice(SquareLattice*, int, int, int);
void initiateNeighbours(Neighbours*, int, int);
void clearSquareLattice(SquareLattice*);
int isValidInteraction(SquareLattice*, int, int);
void labelCluster(int*, int*, SquareLattice*);
int unionFind(int*, int, int);
int getInterfacialLengthWith(int*, SquareLattice*, int, int);




void initiateSquareLattice(SquareLattice *l, int dimension, int L, int B){
    int size = pow(L, dimension);
    l -> dimension = dimension;
    l -> size = size;
    l -> neighbours = safeMAlloc(sizeof(Neighbours));
    initiateNeighbours(l->neighbours, dimension, L);
    l -> boundryCondition = B; 
}

void initiateNeighbours(Neighbours *neighbours, int dimension, int L){

    neighbours -> size = pow(L, dimension);
    neighbours -> neighboursCount = pow(2, dimension);
    neighbours -> matrix = safeMAlloc(neighbours -> size * sizeof(int*));
    
    for(int i = 0; i < neighbours->size; i++){
        neighbours -> matrix[i] = safeMAlloc(neighbours->neighboursCount * sizeof(int));
    }
    for(int i = 0; i < neighbours->size; i++){
        neighbours -> matrix[i][0] = i % L == (L - 1) ? i - L + 1 : i + 1;
        neighbours -> matrix[i][1] = i % L == 0       ? i + L - 1 : i - 1;
        if(dimension > 1){
            int N = pow(L, 2);
            neighbours -> matrix[i][2] = i < L        ? i + N - L : i - L;
            neighbours -> matrix[i][3] = i >= (N - L) ? i - N + L : i + L;  
        }
    }
    
}

void clearSquareLattice(SquareLattice *lattice){
    if(lattice == NULL || lattice->neighbours == NULL) return;
    for(int i = 0; i < lattice->neighbours->size; i++){
        free(lattice->neighbours->matrix[i]);
    }
    free(lattice->neighbours->matrix);
    free(lattice->neighbours);
    lattice->neighbours = NULL;
}


/************isValidInteraction()*************
    Verify if the interaction is valid
    based on the boundry condition
    dim = 1:
        0: ring (periodic);
        1: stripe (non-periodic)
        
    dim = 2:
        0: torus (periodic)
	    1: vertical cilindric (periodic horizontally)
		2: horizontal cilindric (periodic vertically)
	    3: square (non-periodic)
	dim = 3:
	    0: periodic	 
*********************************************/
int isValidInteraction(SquareLattice* lattice, int s1, int s2){
    if(lattice->boundryCondition == 0) return 1;
    if(lattice->dimension >= 3) return 1;
    int L = lattice->L;
    if(lattice->dimension == 1){
        return !(abs(s1 - s2) > 1); // on stripes neighbours have distance equals to 1 always
    } else if(lattice->dimension == 2){
        int x1 = s1 % L, x2 = s2 % L;
        int y1 = s1 / L, y2 = s2 / L;
        int f1 = abs(y1 - y2) <= 1;
        int f2 = abs(x1 - x2) <= 1; // same line or non periodic neighbours
        
        if(lattice->boundryCondition == 1){
            return f1;
        } else if(lattice->boundryCondition == 2){
            return f2;
        } else {
            return f1 && f2;
        }
    }
    return 0;
}



/***************labelCluster()****************
    Creates the Hoshen-Kopelman label list
    label[i] returns the index of the cluster
    the index of the cluster is the same as
    the index of the last agent on that cluster
    
*********************************************/
void labelCluster(int *label, int *s, SquareLattice* lattice){
    for(int i = 0; i < lattice->size; i++){
        label[i] = i;    
    }
    int L = lattice->L;
    for(int i = 0; i < lattice->size; i++){
        int ri = lattice->neighbours->matrix[i][0];
        int di = lattice->neighbours->matrix[i][3];
        if(s[i] == s[ri]){// CHECKING WITH THE RIGHT
            if(i % L != L - 1 || isValidInteraction(lattice, s[i], s[ri])){ 
            label[i] = unionFind(label, i, ri);
            }
        }
        if(s[i] == s[di]){// CHECKING WITH THE BOTTOM
            if(i % L != L - 1 || isValidInteraction(lattice, s[i], s[di])){ 
            label[i] = unionFind(label, i, di);
            }
        }

    }
    return;
}


/***************labelCluster()****************
    Returns the last index of the cluster;
*********************************************/
int unionFind(int *lab, int i1_0, int i2_0){
    int i1 = i1_0, i2 = i2_0, I;
    while(i1 != lab[i1]){
        i1 = lab[i1];    
    }
    while(i2 != lab[i2]){
        i2 = lab[i2];
    }
    I = i1;
    if(i1 < i2) I = i2;
    for(int i = 0; i <= I; i++){
    	if(lab[i] == i1 || lab[i] == i2){
    		lab[i] = I;
    	}
    }
    return I;
}


/**********getInterfacialLengthWith()*********
	Measure the length of the cluster(l1) with
	the other cluster(l2), ONLY THE INTERFACE
	WITH l1 and l2. IF THERES ANOTHER CLUSTER
	WITH INTERFACE WITH L1 
	Uses the H-K algorithm. Sen
	sitive to boundry conditions. 
*********************************************/

int getInterfacialLengthWith(int* lab, SquareLattice* lattice, int l1, int l2){
	int li, lj, length = 0;
	for(int i = 0; i < lattice->size; i++){
		li = lab[i];
		if(li != l1) continue;
		for(int j = 0; j < lattice->neighbours->neighboursCount; j++){
			if(!isValidInteraction(lattice, i, lattice->neighbours->matrix[i][j])) continue;
			lj = lab[lattice->neighbours->matrix[i][j]];
			if(lj != l2) continue;
			length += 1;
		}
	}
	return length;
}





