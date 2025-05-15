#include <lplattice.h>
#include <lpfhelper.h>
#include <ploting.h>
#include <string.h>

void setup(void);
void mergeOldValues(void);
void memoryAllocation(void);
void buildInitialConfiguration(void);
void buildInitialMobileMatrix(void);
void initiate(void);
void take1DMeasures(double);
void take2DMeasures(double);
void takeMeasures(double);
void clear(void);
void singleInteraction(void);
void updateMobileMatrix(int);
int isMobile(int);
void locMatrixUpdate(int, int, int*, int*);
void openFile(void);
void temporalEvolution(void);
void reinforcementEvolution(void);
void lengthEvolution(void);
float correlationTime(float);
float timeForWStat(float);
void evolveSystemTo(double);
void printAll(void);
void writeInstructions(FILE*, int);
void write2DInstructions(FILE*, int);
void write1DInstructions(FILE*, int);

#define SEED                0 // 0: random

//==================Lattice Config.==================//
#define LSIZE               128
#define DIM                 2

//====================Generic.Evolution================//
#define DETA                1e-2

//==================Temp. Evolution==================//
#define TMAX                1e6
#define TMIN                1e0  // cant be 0 for scale log

//==================Length. Evolution================//
#define LMAX                64	
#define LMIN                8 

//==================DEta. Evolution==================//
#define DETAMIN             1e-4
#define DETAMAX             1e0

//===============Initial Configuration===============//
#define C                   2           // 0: Random
                                        // 1: Half(Z/Z)
                                        // 2: Circular


//===============Measures Configuration==============//

#define MEASURES            60          // Number of evolution tics
#define LOOPS               1e1         // Number of loops(only valid for non temporal evolutions)
#define SCALE               1           // 0: Linear
                                        // 1: LOG
#define MODE                1           // 0: Number of zealots
                                        // 1: cluster
#define EVOLUTION           0           // 0: Temporal
                                        // 2: Length Reinforcement
                                        // 2: Spatial

//=================Boundry Condition=================//
#define FIXED_BORDERS       0

#define B				    0			    //"b" = "boundry condition"
        									//0 -> periodic
											//1 -> dobrushin vertical(periodic horizontally)
											//2 -> dobrushin horizontal(periodic vertically)
											//3 -> square
#define PLOT                0


int N = pow(LSIZE, DIM);
SquareLattice lattice;
unsigned int seed;
FILE *fp1, *fCompl;

double *eta, dEta;
double timeT;
// always use s and z to measures;
int *s, *post_s, *z, *post_z;

//===============DIM: 2; MODE: 1=====================//
//int nClS_M, lengthCl1S_M, lengthCl2S_M, nClSZ_M, lengthCl1SZ_M, lengthCl2SZ_M;

// mobile its a matrix of mobile agents
// mobileCount it its size
// locMobile its a matrix of mobile matrix agents index
// locMobile[i] = index of agent i in mobile matrix
int mobileCount, *mobile, *locMobile;
float genMeasure;

//

int main(){
    initiateSquareLattice(&lattice, DIM, LSIZE, B);
    setup();
    initiate();
    clear();
    return 0;
}




void setup(void){
    memoryAllocation();
}


void mergeOldValues(void){
    for(int i = 0; i < N; i++){
        s[i] = post_s[i];
        z[i] = post_z[i];
    }

    return;
}

void memoryAllocation(void){
    s = safeMAlloc(N * sizeof(int));
    post_s = safeMAlloc(N * sizeof(int));
    z = safeMAlloc(N * sizeof(int));
    post_z = safeMAlloc(N * sizeof(int));
    mobile = safeMAlloc(N * sizeof(int));
    locMobile = safeMAlloc(N * sizeof(int));
    eta = safeMAlloc(N * sizeof(double));

}

void buildInitialConfiguration(void){
    switch(C){
        case 0:
            for(int i = 0; i < N; i++){
                post_s[i] = randomInt(2);
                post_z[i] = 1;
                eta[i] = 1.0;
            }
            break;
        case 1:
            for(int i = 0; i < N; i++){
                post_s[i] = i < (N / 2) ? 0 : 1;
                post_z[i] = 1;
                eta[i] = 1.0;
            }
            break;
        case 2:
        	float x0 = LSIZE / 2.0, y0 = x0;
            for(int i = 0; i < N; i++){
                int x = i % LSIZE, y = i / LSIZE;
                float distance = sqrt(pow(x - x0, 2) + pow(y - y0, 2));

                post_s[i] = distance <= (float) lattice.L / 4.0 ? 1 : 0;
                post_z[i] = 1;
                eta[i] = 1.0;
            }
            break;

    }
    mergeOldValues();
    if(PLOT) printAll();
    buildInitialMobileMatrix();

}
void buildInitialMobileMatrix(void){
    mobileCount = 0;
    for(int i = 0; i < N; i++){
        if(isMobile(i)){
            mobile[mobileCount] = i;
            locMobile[i] = mobileCount;
            mobileCount++;
        } else {
            mobile[N - 1 + (mobileCount - i)] = i;
            locMobile[i] = N - 1 + (mobileCount - i);
        }
    }
}
void initiate(void){
    switch(EVOLUTION){
        case 0: //temporalEvolution
            temporalEvolution();
            fclose(fp1);
            break;
        case 1:
            reinforcementEvolution();
            break;
        case 2:
            lengthEvolution();
            break;
    }
}

void printAll(){
    for(int i = 0; i < N; i++){
        printf("%d ", s[i] + 2 * z[i]);
        if(DIM == 2 && i % LSIZE == LSIZE - 1) printf("\n");
    }

    printf("\n\n");
}



void temporalEvolution(void){
    openFile();
	float *time_arr = safeMAlloc(MEASURES * sizeof(float));
	if(SCALE){
		geomProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
	} else {
		linearProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
	}


	timeT = 0.0;
	double nextTime = 0.0;
    buildInitialConfiguration();
    
    dEta = (double)DETA;
    for(int i = 0; i < MEASURES; i++){
    	nextTime = (double)time_arr[i];
        evolveSystemTo(nextTime);
    	takeMeasures(nextTime);
    	//printAll();
    }
    free(time_arr);
}

void lengthEvolution(void){
	float *lengthArr = safeMAlloc((int)MEASURES * sizeof(float));
	if(SCALE){
		geomProgression(lengthArr, (float)LMIN, (float)LMAX, (int)MEASURES);
	} else {
		linearProgression(lengthArr, (float)LMIN, (float)LMAX, (int)MEASURES);
  	}
  	fCompl = safeOpen("data_MEAN", ".dat");
	writeInstructions(fCompl, 1);
  	dEta = (double)DETA;
  	for(int i = 0; i < (int)MEASURES; i++){
    	
  		clear();
		int length = (int)lengthArr[i];
		
		initiateSquareLattice(&lattice, DIM, length, B);
		N = (int)pow(length, DIM);
  		openFile();
  		//printf("S %d %d %d\n", length, lattice.L, N);
		setup();

		timeT = 0;
		buildInitialConfiguration();
		
		double nextTime = timeForWStat(dEta);
		evolveSystemTo(nextTime);
		for(int l = 0; l < (int)LOOPS; l++){
		    nextTime += correlationTime(dEta);
		    evolveSystemTo(nextTime);
		    takeMeasures(nextTime);
		}
		genMeasure = 0;
		fclose(fp1);
		fprintf(fCompl, "%.5f %.5f\n", dEta, genMeasure / LOOPS);
		fflush(fCompl);
  	}
  	free(lengthArr);
}

void reinforcementEvolution(void){
    float *deltaEtaArr = safeMAlloc((int)MEASURES * sizeof(float));
    if(SCALE){
		    geomProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, (int)MEASURES);
    } else {
		    linearProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, (int)MEASURES);
	  }


	fCompl = safeOpen("data_MEAN", ".dat");
	writeInstructions(fCompl, 1);
    for(int i = 0; i < (int)MEASURES; i++){
        openFile();
        genMeasure = 0;
        timeT = 0;
        dEta = deltaEtaArr[i];
        buildInitialConfiguration();

        double nextTime = timeForWStat(dEta);
        evolveSystemTo(nextTime);

        for(int l = 0; l < (int)LOOPS; l++){
            nextTime += correlationTime(dEta);
            evolveSystemTo(nextTime);
            takeMeasures(nextTime);
        }
        fclose(fp1);
        fprintf(fCompl, "%.5f %.5f\n", dEta, genMeasure / LOOPS);
        fflush(fCompl);
    }
    fclose(fCompl);
    free(deltaEtaArr);
}

float timeForWStat(float dEta){
    return 10 / dEta;
}
float correlationTime(float dEta){
    return 10 / dEta;
}
void evolveSystemTo(double nextTime){
    while(timeT < nextTime){
        if(mobileCount == 0){
            mergeOldValues();
            break;
        }
        timeT += 1. / (double)mobileCount;
    	if(timeT >= nextTime){
    	    mergeOldValues();
    	}
    	singleInteraction();
    }
}


void take1DMeasures(double nextTime){
    switch(MODE){
        case 0:{
            int num_z = 0;
            int num_s1 = 0;
            int num_z1 = 0;
            int lastZ0 = 0, foundS1 = 0;
            for(int i = 0; i < N; i++){
                num_s1 += s[i];
                num_z1 += s[i] * z[i];
                num_z  += z[i];
            	if(s[i] == s[0]){ //
            	    if(z[i] && !foundS1){
            	        lastZ0 = i;
            	    }
            	} else {
            	    foundS1 = 1;
            	}

            }
            fprintf(fp1,"%.2f %d %d %d %d\n", nextTime, num_s1, num_z1, num_z, lastZ0);

            break;
        }
        case 1:{
            int nz = 0;
            for(int i = 0; i < N; i++){
                nz += z[i];

            }
            genMeasure += nz;
            fprintf(fp1,"%.2f %d\n", nextTime, nz);
            break;


        }
    }

}
void take2DMeasures(double nextTime){
	if(N <= 1e2 && PLOT) printAll(); 
    switch(MODE){
		case 0:
            int num_z = 0, num_s1 = 0, num_z1 = 0;
            for(int i = 0; i < N; i++){
                num_s1 += s[i];
                num_z  += z[i];
                num_z1 += s[i] * z[i];
            }
            fprintf(fp1,"%.2f %d %d %d\n", nextTime, num_s1, num_z1, num_z);
            break;
        case 1: //clusters length
        	
			int *matrix_SZ = safeMAlloc(N * sizeof(int));
         	for(int i = 0; i < N; i++){
            	matrix_SZ[i] = s[i] + 2 * z[i];
          	}

	        int labCl1SZ, labCl2SZ, nClSZ = 0;//biggestClustersLabels and number of Clusters of S
    	    int *labSZ = safeMAlloc(N * sizeof(int)); //H-K label matrix of S;
        	int *sizeSZ = safeMAlloc(N * sizeof(int)); // Sizes of clusters;
		    labelCluster(labSZ, matrix_SZ, &lattice);
		    clusterSizeInfo(labSZ, sizeSZ, &nClSZ, &labCl1SZ, &labCl2SZ, &lattice);
			int lengthCl1SZ = getInterfacialLength(labSZ, &lattice, labCl1SZ); //length of biggest cluster of S
            int lengthCl2SZ = getInterfacialLength(labSZ, &lattice, labCl2SZ); //length of SECOND biggest cluster of S


          	int labCl1S, labCl2S, nClS = 0;//biggestClustersLabels and number of Clusters of S
          	int *labS = safeMAlloc(N * sizeof(int)); //H-K label matrix of S;
          	int *sizeS = safeMAlloc(N * sizeof(int)); // Sizes of clusters;
          	labelCluster(labS, s, &lattice);
          	clusterSizeInfo(labS, sizeS, &nClS, &labCl1S, &labCl2S, &lattice);
          	int lengthCl1S = getInterfacialLength(labS, &lattice, labCl1S); //length of biggest cluster of S
          	int lengthCl2S = getInterfacialLength(labS, &lattice, labCl2S); //length of SECOND biggest cluster of S

			
          fprintf(fp1, "%.2f %d %d %d %d %d %d\n", nextTime, nClS, lengthCl1S, lengthCl2S, nClSZ, lengthCl1SZ, lengthCl2SZ);
		
          free(labS);  free(matrix_SZ); free(sizeS);
          free(labSZ); free(sizeSZ);
    }
}

void takeMeasures(double nextTime){
    if(DIM == 1){
        take1DMeasures(nextTime);
    } else {
        take2DMeasures(nextTime);
    }
    fflush(fp1);
}


void clear(void){
    free(s); free(post_s);
    free(z); free(post_z);
    free(mobile); free(locMobile);
    clearSquareLattice(&lattice);
}

//==================Virtual Environment==================//
//      It acts wit post_s and post_z, virtual versions
//      of s and z matrixes. When the time of the simu-
//      lation gets closer to the time for the measure
//      the virtual values from post_s and post_z are
//      merged for s and z. Never use post_s and post_z
//      to take measures, values from this matrixes do not
//      correspond for the values of timeT.
void singleInteraction(void){
    int changedState = 0;


    int pPos = mobile[randomInt(mobileCount)];
    int vPos = lattice.neighbours->matrix[pPos][randomInt(lattice.neighbours->neighboursCount)];

	if(!isValidInteraction(&lattice, pPos, vPos)) {
        //printf("%d %d\n", pPos, vPos);
        return;
	}
	if(FIXED_BORDERS) {
	    int on1DBorder = (pPos == 0 || pPos == LSIZE - 1);
	    if(DIM == 1 && on1DBorder) {
	        if(C == 2) {
	            if(pPos == 0)return;
	        } else {
	            return;
	        }

	    } else if(DIM == 2){
	        int on2DBorder = pPos < LSIZE || pPos >= N - LSIZE;
	        on2DBorder |= (B == 3 && (pPos % LSIZE == 0 || pPos % LSIZE == LSIZE - 1));
	        if(on2DBorder) {
	            //printf("%d\n", pPos);
	            return;
	        }
	    }

	}
    int reinforcementInteraction = post_s[pPos] == post_s[vPos];

    if(reinforcementInteraction){
        eta[pPos] += dEta;
        if(eta[pPos] >= 1.){
            post_z[pPos] = 1;
            changedState = 1;
        }
    } else {
        if(post_z[pPos] == 1){
            post_z[pPos] = 0;
            eta[pPos] = 0.0;
        } else {
            post_s[pPos] = post_s[vPos];
            eta[pPos] = 0.0;
            changedState = 1;
        }
    }

    if(changedState){
        updateMobileMatrix(pPos);
        for(int i = 0; i < lattice.neighbours->neighboursCount; i++){
            updateMobileMatrix(lattice.neighbours->matrix[pPos][i]);
        }
    }
    return;
}
void updateMobileMatrix(int i){

    if(isMobile(i)){

        if(locMobile[i] >= mobileCount){ // is mobile but its not on the mobile list

            // change pos of i with the last element of the mobile

            locMatrixUpdate(i, mobile[mobileCount], mobile, locMobile);
            mobileCount++;
        }
    } else {
        if(locMobile[i] < mobileCount){ // its not mobile but is on the llist
            // change pos of i with the last element of the mobile
            locMatrixUpdate(i, mobile[mobileCount - 1], mobile, locMobile);

            mobileCount--;
        };
    }
}
int isMobile(int i){
    if(post_z[i] == 1){
        for(int j = 0; j < lattice.neighbours->neighboursCount; j++){
            int neigIndex = lattice.neighbours->matrix[i][j];
            if(post_s[neigIndex] != post_s[i]){
                return 1;
            }
        }
        return 0;
    }
    return 1;
}
void locMatrixUpdate(int s1, int s2, int *mat, int *locMat){
    int tempS;

    mat[locMat[s1]] = s2;
    mat[locMat[s2]] = s1;

    tempS = locMat[s1];
    locMat[s1] = locMat[s2];
    locMat[s2] = tempS;
    return;


}

void openFile(){
    seed = setupRandom(SEED);
    char name[100];
    switch(MODE){
      	case 0:
      		sprintf(name, "data_pvm_Ev%d_C%d_B_%d_L%d", MODE, C, B, LSIZE);
      		break;
      	case 1:
      		sprintf(name, "data_pvm_Ev%d_C%d_B_%d_L%d", MODE, C, B, LSIZE);
            break;
      	case 2:
      		sprintf(name, "data_pvm_Ev%d_C%d_B_%d_dE%.4f", MODE, C, B, dEta);
      		break;
    }
    fp1 = safeSeedOpen(name, ".dat", &seed, SEED != 0);
    writeInstructions(fp1, 0);
    fflush(fp1);

}

void writeInstructions(FILE* f, int isCompl){

	fprintf(f, "# Persistent Voter Model: pvm.c\n");
	fprintf(f, "# Data generated by: L. Possamai\n");
	fprintf(f, "# Seed: %d\n", seed);
	fprintf(f, "# Dimension: %d\n", DIM);
	if(SCALE){
	    fprintf(f, "# Log Measures: %d\n", (int)MEASURES);
	} else {
		fprintf(f, "# Linear Measures: %d\n", (int)MEASURES);
	}
    switch(EVOLUTION){
    	case 0:
    		fprintf(f, "# dEta = %.5f\n", DETA);
			fprintf(f, "# L = %d\n", (int)LSIZE);
			fprintf(f, "# t: [%d, %d]\n", (int) TMIN, (int)TMAX);
    		break;
    	case 1:
			fprintf(f, "# L = %d\n", (int)LSIZE);
			if(isCompl){
			    fprintf(f, "# loops: %d\n", (int)LOOPS);
				fprintf(f, "# dEta: [%.5f, %.5f]\n", DETAMIN, DETAMAX);
			} else {
				fprintf(f, "# dEta: %.5f\n", (double)dEta);
			}
			break;
		case 2:
    		fprintf(f, "# dEta = %.5f\n", DETA);
			if(isCompl){
			    fprintf(f, "# loops: %d\n", (int)LOOPS);
				fprintf(f, "# L: [%d, %d]\n", (int)LMIN, (int)LMAX);
			} else {
				fprintf(f, "# L: %d\n", (int)lattice.L);
			}
			break;
    }
    if(DIM == 1){
        write1DInstructions(f, isCompl);
    } else if(DIM == 2){
        write2DInstructions(f, isCompl);
    }
}

void write2DInstructions(FILE* f, int isCompl){
    switch(MODE){
        case 0:
			fprintf(f, "#  t nS0 nZ0 nZ\n");
			break;    
		case 1:
			fprintf(f, "#  t nCl_S 1S_len 2S_len nCl_SZ 1SZ_len 2SZ_len\n");
			break;
    }
}
void write1DInstructions(FILE* f, int isCompl){

	switch(MODE){
		case 0:
			fprintf(f, "#  t N_S0 N_Z0 N_Z lastZ0\n");
			break;
		case 1:
		    if(isCompl){
			    fprintf(f, "#  t <nZ>\n");
			} else {
			    fprintf(f, "#  t nZ\n");
			}
			break;
	}
}
