//#include <monte_carlo.h>
#include <lattice.h>
#include <fhelper.h>
#include <ploting.h>
#include <string.h>

void setup(void);
void mergeOldValues(void);
void memoryAllocation(void);
void buildInitialConfiguration(void);
void buildInitialMobileMatrix(void);
void initiate(void);
void takeMeasures(double);
void clear(void);
void singleInteraction(void);
void updateMobileMatrix(int);
int isMobile(int);
void locMatrixUpdate(int, int, int*, int*);
void openFile(void);
void writeInstructions(FILE*, int);
void temporalEvolution(void);
void reinforcementEvolution(void);
float correlationTime(float);
float timeForWStat(float);
void evolveSystemTo(float);
void printAll(void);

#define SEED                0 // 0: random
#define L                   150
#define DIM                 1

//==================Temp. Evolution==================//
#define DETA                1e0
#define TMAX                1e6
#define TMIN                0  // cant be 0 for scale log
//==================DEta. Evolution==================//
#define DETAMIN             1e-4
#define DETAMAX             1e0
#define LOOPS               1e0
//===================================================//

#define MEASURES            2

#define C                   1           // 0: Random
                                        // 1: Half
#define SCALE               1           // 0: Linear
                                        // 1: LOG
#define MODE                1           // 0: Temporal
#define FIXED_BORDERS       1

#define B					0			    //"b" = "boundry condition"  
        									//0 -> periodic
											//1 -> dobrushin vertical(periodic horizontally)
											//2 -> dobrushin horizontal(periodic vertically)
											//3 -> square
#define PLOT                0


int N = pow(L, DIM);
SquareLattice lattice;
unsigned int seed;
FILE *fp1, *fCompl;

double *eta, dEta;
double timeT; 
// always use post_s and post_z to measures;
int *s, *post_s, *z, *post_z;

// mobile its a matrix of mobile agents
// mobileCount it its size
// locMobile its a matrix of mobile matrix agents index
// locMobile[i] = index of agent i in mobile matrix
int mobileCount, *mobile, *locMobile;
float genMeasure;

int main(){
    initiateSquareLattice(&lattice, DIM, L);
    if(PLOT){
        startLatticeGif(L);
        char **l;
		l = smalloc(4 * sizeof(char*));
		l[0] = "0xEDC001";
		l[1] = "0xFFE135";
		l[2] = "0x4B0082"; 
		l[3] = "0x663399";
		setColorPallete(l, 4);
    }
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
    s = smalloc(N * sizeof(int));
    post_s = smalloc(N * sizeof(int));
    z = smalloc (N * sizeof(int));
    post_z = smalloc(N * sizeof(int));
    mobile = smalloc(N * sizeof(int));
    locMobile = smalloc(N * sizeof(int));
    eta = smalloc(N * sizeof(double));
    
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
    }
    mergeOldValues();
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
    switch(MODE){
        case 0: //temporalEvolution
            temporalEvolution();
            fclose(fp1);
            break;
        case 1:
            reinforcementEvolution();
            break;
    }
}

void printAll(){
    printf("\n%d mobiles: ", mobileCount);
    for(int i = 0; i < mobileCount; i++){
        printf("%d ", mobile[i]);
    }
    printf("\n");
    for(int i = 0; i < N; i++){
        printf("%d ", post_s[i] + 2 * post_z[i]);
    }
    printf("\n");
    for(int i = 0; i < N; i++){
        printf("%d ", isMobile(i));
    }

    printf("\n");
}

void temporalEvolution(void){
    openFile();
	float *time_arr = smalloc(MEASURES * sizeof(float));
	if(SCALE){
		geomProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
	} else {
		linearProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
	}
	
	
	timeT = 0.0;
	double nextTime = 0.0;
    buildInitialConfiguration();
    dEta = DETA;
    for(int i = 0; i < MEASURES; i++){
    	nextTime = time_arr[i];
        evolveSystemTo(nextTime);
    	takeMeasures(nextTime);
    	
    	if(PLOT){
    	    int* s0 = smalloc(N * sizeof(int));
    	    for(int i = 0; i < N; i++){
    	        s0[i] = (int)(2 * post_s[i]) + post_z[i];
    	    }
    	    printLine(s0, L);
    	    char timeLabel[100];
    	    sprintf(timeLabel, "%.2f", nextTime);
    	    printLabelAt(" t = ", timeLabel, 0.5, 0.95);
    	    
    	    free(s0);
    	}
    }
    free(time_arr);
}
void reinforcementEvolution(void){
    float *deltaEtaArr = smalloc(MEASURES * sizeof(float));
    if(SCALE){
		geomProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, MEASURES);
	} else {
		linearProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, MEASURES);
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
void evolveSystemTo(float nextTime){
    while(timeT < nextTime && mobileCount != 0){
        mergeOldValues();
        timeT += 1. / mobileCount;
    	singleInteraction();
    }
}


void takeMeasures(double nextTime){
    switch(MODE){
        case 0:{
            int nz = 0;
            int s1 = 0;
            int z1 = 0;
            
            int lastZ0 = 0, firstZ1 = L - 1;
            for(int i = 0; i < N; i++){
            	if(s[i] == 1){
            		z1 += z[i];
            		s1++;
            		
                }
                if(z[i]){
                    if(s[i] == s[0]){ // last zealot with same opinion as i = 0
                        lastZ0 = i;
                    } else if(i < firstZ1) { // first zealot with different opinion of i=0;
                        firstZ1 = i;
                    }
                }
                nz += z[i];
            }
            //clusterNumber(&ncl,&per);
            //verificarCruzamentoInterface();
            fprintf(fp1,"%.2f  %d %d %d %d %d\n", nextTime, s1, z1, nz, N - nz, firstZ1 - lastZ0 - 1);
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

void clear(void){
    free(s); free(post_s);
    free(z); free(post_z);
    free(mobile); free(locMobile);
}

void singleInteraction(void){
    int changedState = 0;

    
    int pPos = mobile[randomInt(mobileCount)];
    int vPos = lattice.neighbours.matrix[pPos][randomInt(lattice.neighbours.neighboursCount)];
    
	if(!validInteraction(&lattice, pPos, vPos, B)) return;
	if(FIXED_BORDERS && (pPos == 0 || pPos == L-1)) return;
    int reinforcementInteraction = s[pPos] == s[vPos];
    
    if(reinforcementInteraction){
        eta[pPos] += dEta;
        if(eta[pPos] >= 1.){
            post_z[pPos] = 1;
            changedState = 1;
        }
    } else {
        if(z[pPos] == 1){
            post_z[pPos] = 0;
            eta[pPos] = 0.0;
        } else {
            post_s[pPos] = s[vPos];
            eta[pPos] = 0.0;
            changedState = 1;
        }
    }

    if(changedState){
        updateMobileMatrix(pPos);
        for(int i = 0; i < lattice.neighbours.neighboursCount; i++){
            updateMobileMatrix(lattice.neighbours.matrix[pPos][i]);
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
        for(int j = 0; j < lattice.neighbours.neighboursCount; j++){
            int neigIndex = lattice.neighbours.matrix[i][j];
            if(post_s[neigIndex] != s[i]){
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
      		sprintf(name, "data_pvm_Ev%d_C%d_B_%d_L%d", MODE, C, B, L);
      		break;
      	case 1:
      		sprintf(name, "data_pvm_Ev%d_C%d_B_%d_L%d", MODE, C, B, L);
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
	
	switch(MODE){	
		case 0:
			fprintf(f, "# dEta: %.5f\n", (double)DETA);
			fprintf(f, "# L: %d\n", L);
			fprintf(f, "#  t  nS1  nZ1  nZ nNormal 1dW\n");
			break;
		case 1:
		
			fprintf(f, "# L: %d\n", L);
		    if(isCompl){
			    fprintf(f, "# dEta: [%.5f, %.5f]\n", DETAMIN, DETAMAX);
			    fprintf(f, "# loops: %d\n", (int)LOOPS);
			    fprintf(f, "#  t <nZ>\n");
			} else {    
			    fprintf(f, "# dEta: %.5f\n", (double)DETA);
			    fprintf(f, "#  t nZ\n");
			
			}
			break;
		case 2:
			//fprintf(fp1, "# dEta = %.5f\n", DELTAETA_INICIAL);
			//fprintf(fp1, "# Tau_e = %d\n", (int)TEMPO_MAX);
			//fprintf(fp1, "# L l_i l_arr d_f\n");
			break;
		case 3:
			//fprintf(fp1, "# dEta = %.5f\n", DELTAETA_INICIAL);
			//fprintf(fp1, "# L = %d\n", LSIZE);
			//fprintf(fp1, "# Tau_e = %d\n", (int)TEMPO_MAX);
			//fprintf(fp1, "#  t int_z int_ar\n");
			break;
	}
}






