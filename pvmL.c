#include <lpfhelper.h>
#include <lplattice.h>
#include <lputil.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
void temporalEvolution(void);
void reinforcementEvolution(void);
void lengthEvolution(void);
void inversionTime1D(void);
double correlationTime(double);
double timeForWStat(double);
void evolveSystemTo(double);
void printAll(void);
void openFile(FILE**, char*);
void w1DBorders(int*, int*);
void writeInstructions(FILE*, int, int);
void write2DInstructions(FILE*, int, int);
void write1DInstructions(FILE*, int, int);
void plotSMatrix(void);
void onBorderWalkEvolution(void);
double etaOfBiased(int, int);
void exportBaseSpinMatrix(void);
void attainConsensus(void);
void initiateGenerics(int);
void printGenerics(FILE*, int);
int mc_biasedwalk(int qual, int end, int* lab, double*);

bool isOnBorder(int);

//==================Lattice Config.==================//

#ifndef EXPORT_MATRIX
#define EXPORT_MATRIX 1
#endif
#ifndef SEED
#define SEED 0 // 123456789 //: random
#endif
//==================Lattice Config.==================//

#define LSIZE 128
#define HSIZE 128
#define DIM 2

//==================Generic.Evolution================//

#define DETA 1E-2

//==================Temp. Evolution==================//

#define TMIN 1e0 // cant be 0 for scale log
#define TMAX 1e4
//==================Length Evolution=================//

#define LMAX 256
#define LMIN 32

//==================DEta. Evolution==================//

#define DETAMIN 1e-2
#define DETAMAX 1e0

//===============Initial Configuration===============//

/*
 * 0: Random
 * 1: Half(Z/Z)
 * 2: Circular
 */
#define C 1

//===============Measures Configuration==============//

// Number of evolution tics
#define MEASURES 30

// Number of loops(only valid for non temporal evolutions)
#define LOOPS 1e3

// 0: Linear
// 1: LOG
#define SCALE 1

// 0: Number of zealots
// 1: Interface Length
// 2: Instantaneous AR Width
#define MODE 0

#define ALL_ZEALOTS 1
/*
 * 0: temporalEvolution
 * 1: reinforcementEvolution
 * 2: lengthEvolution
 * 3: inversionTime1D
 * 4: onBorderWalkEvolution
 */
#define EVOLUTION 0

//=================Boundry Condition=================//
#define FIXED_BORDERS 0

/*
 * "b" = "boundry condition"
 * 0 -> periodic
 * 1 -> dobrushin vertical(periodic horizontally)
 * 2 -> dobrushin horizontal(periodic vertically)
 * 3 -> square
 */
#define B 1

// ==================================================//

int N;
Dimension dimension;
Lattice lattice;
unsigned int seed = -1;
FILE *fp1, *fCompl, *fMatrix;

double *eta, *post_eta, dEta;
double timeT;
// always use s and z to measures;
int *s, *post_s, *z, *post_z;
//===============DIM: 2; MODE: 1=====================//
// int nClS_M, lengthCl1S_M, lengthCl2S_M, nClSZ_M, lengthCl1SZ_M, lengthCl2SZ_M;

// mobile its a matrix of mobile agents
// mobileCount it its size
// locMobile its a matrix of mobile matrix agents index
// locMobile[i] = index of agent i in mobile matrix
int mobileCount, *mobile, *locMobile;

int nZ = 0, nS1 = 0, nZ1 = 0, nClean = 0;
double gen0 = 0, gen1 = 0, gen2 = 0, gen3 = 0;
double* generics = NULL;
int genericsLength = 0;

int main() {
    createDimension(&dimension, DIM, (size_t[]){LSIZE, HSIZE});
    N = dimension.size;
    createLattice(&lattice, &dimension, B);

    if (EXPORT_MATRIX) {
        fMatrix = fopen("matrix.dat", "wb");
    }
    setup();
    initiate();
    clear();
    if (EXPORT_MATRIX && fMatrix != NULL) {
        fclose(fMatrix);
    }
    return 0;
}

void initiate(void) {

    switch (EVOLUTION) {
        case 0: // temporalEvolution
            temporalEvolution();
            break;
        case 1:
            reinforcementEvolution();
            break;
        case 2:
            lengthEvolution();
            break;
        case 3:
            //            inversionTime1D();
            break;
        case 4:
            //            onBorderWalkEvolution();
            break;
        case 5:
            //            attainConsensus();
            break;
    }
}

void setup(void) {
    memoryAllocation();
}

void mergeOldValues(void) {
    for (int i = 0; i < N; i++) {
        s[i] = post_s[i];
        z[i] = post_z[i];
        eta[i] = (double)post_eta[i];
    }

    return;
}

void memoryAllocation(void) {
    s = safeMAlloc(N * sizeof(int));
    post_s = safeMAlloc(N * sizeof(int));
    z = safeMAlloc(N * sizeof(int));
    post_z = safeMAlloc(N * sizeof(int));
    mobile = safeMAlloc(N * sizeof(int));
    locMobile = safeMAlloc(N * sizeof(int));
    eta = safeMAlloc(N * sizeof(double));
    post_eta = safeMAlloc(N * sizeof(double));
}

bool isOnBorder(int i) {

    bool f = false;
    for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
        // IS A BORDER IF A INTERACTION IS BLOCKED
        if (!isValidInteraction(&lattice, i, lattice.neighbours->matrix[i][j])) {
            f = true;
            break;
        }
    }
    return f;
}

void buildInitialConfiguration(void) {
    switch (C) {
        case 0:
            for (int i = 0; i < N; i++) {
                post_s[i] = randomInt(2);
                post_z[i] = ALL_ZEALOTS;
                post_eta[i] = ALL_ZEALOTS;
                if (FIXED_BORDERS && isOnBorder(i)) {
                    post_z[i] = 1;
                    post_eta[i] = 1.0;
                }
            }
            break;
        case 1:
            for (int i = 0; i < N; i++) {
                post_s[i] = i < (N / 2) ? 0 : 1;
                post_z[i] = ALL_ZEALOTS;
                post_eta[i] = ALL_ZEALOTS;
                if (FIXED_BORDERS && isOnBorder(i)) {
                    post_z[i] = 1;
                    post_eta[i] = 1.0;
                }
            }
            break;
        case 2: {
            if (DIM == 1)
                printf("\n Initial configuration not valid for given dimension \n");
            int i0 = LSIZE % 2 == 0 ? (N + LSIZE) / 2 : N / 2;
            int x0 = i0 % LSIZE, y0 = i0 / LSIZE;
            for (int i = 0; i < N; i++) {
                int x = i % LSIZE, y = i / LSIZE;
                float distance = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
                post_s[i] = distance <= (float)LSIZE / 4.0f ? 0 : 1;
                post_z[i] = ALL_ZEALOTS;
                post_eta[i] = ALL_ZEALOTS;
            }
            break;
        }
    }
    mergeOldValues();
    buildInitialMobileMatrix();
}
void buildInitialMobileMatrix(void) {
    mobileCount = 0;
    for (int i = 0; i < N; i++) {
        if (isMobile(i)) {
            mobile[mobileCount] = i;
            locMobile[i] = mobileCount;
            mobileCount++;
        } else {
            mobile[N - 1 + (mobileCount - i)] = i;
            locMobile[i] = N - 1 + (mobileCount - i);
        }
    }
}

void exportMatrix(const int* s, int sizeN) {
    for (int i = 0; i < sizeN; i++) {
        fprintf(fMatrix, "%d ", s[i]);
    }
    fprintf(fMatrix, "\n");
}

void exportBaseSpinMatrix() {
    int* arr = safeMAlloc((int)N * sizeof(int));

    for (int i = 0; i < N; i++) {
        arr[i] = s[i] + 2 * z[i];
    }
    exportMatrix(arr, N);
    free(arr);
}

void temporalEvolution(void) {
    openFile(&fp1, "temporal");
    writeInstructions(fp1, EVOLUTION, 0);
    float* time_arr = safeMAlloc(MEASURES * sizeof(float));
    if (SCALE) {
        geomProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
    } else {
        linearProgression(time_arr, (float)TMIN, (float)TMAX, MEASURES);
    }

    timeT = 0.0;
    double nextTime = 0.0;
    dEta = (double)DETA;
    buildInitialConfiguration();

    if (EXPORT_MATRIX) {
        exportBaseSpinMatrix();
    }
    for (int i = 0; i < MEASURES; i++) {
        nextTime = (double)time_arr[i];
        evolveSystemTo(nextTime);
        takeMeasures(nextTime);
    }
    free(time_arr);
    fclose(fp1);
}

void attainConsensus(void) {
    openFile(&fp1, "consensus_MEAN");
    writeInstructions(fp1, 0, 0);
    fprintf(fp1, "# dEta tau\n");
    float* dEta_arr = safeMAlloc(MEASURES * sizeof(float));
    if (SCALE) {
        geomProgression(dEta_arr, (float)DETAMIN, (float)DETAMAX, MEASURES);
    } else {
        linearProgression(dEta_arr, (float)DETAMIN, (float)DETAMAX, MEASURES);
    }
    for (int i = 0; i < MEASURES; i++) {
        dEta = (double)dEta_arr[i];
        int l = 0;
        float tau = 0;
        while (l < LOOPS) {

            buildInitialConfiguration();
            timeT = 0.0;
            int mag = sumIntArray(s, N);
            while (mag != 0 && mag != N) {
                evolveSystemTo(timeT + (1.0 / mobileCount));
                mag = sumIntArray(s, N);
            }
            tau += timeT;
            l++;
        }
        fprintf(fp1, "%.5f %.2f\n", dEta, tau / (float)LOOPS);
    }
    fclose(fp1);
    free(dEta_arr);
}

void lengthEvolution(void) {
    openFile(&fCompl, "length_MEAN");
    writeInstructions(fCompl, MODE, 1);
    fflush(fCompl);
    float* lengthArr = safeMAlloc((int)MEASURES * sizeof(float));
    if (SCALE) {
        geomProgression(lengthArr, (float)LMIN, (float)LMAX, (int)MEASURES);
    } else {
        linearProgression(lengthArr, (float)LMIN, (float)LMAX, (int)MEASURES);
    }
    dEta = (double)DETA;
    for (int i = 0; i < (int)MEASURES; i++) {
        clear();
        int length = (int)round(lengthArr[i]);
        createDimension(&dimension, DIM, (size_t[]){LSIZE, HSIZE});
        createLattice(&lattice, &dimension, B);
        N = (int)pow(length, DIM);
        char name[100];
        sprintf(name, "length_L%d", length);
        openFile(&fp1, name);
        writeInstructions(fp1, MODE, 0);
        setup();

        timeT = 0;
        buildInitialConfiguration();

        double nextTime = timeForWStat(dEta);
        evolveSystemTo(nextTime);
        for (int l = 0; l < (int)LOOPS; l++) {
            nextTime += correlationTime(dEta);
            evolveSystemTo(nextTime);
            takeMeasures(nextTime);
        }
        fprintf(fCompl, "%d", length);
        printGenerics(fCompl, 1);

        fclose(fp1);
    }
    fclose(fCompl);
    free(lengthArr);
}

void reinforcementEvolution(void) {
    openFile(&fCompl, "reinforcement_MEAN");
    writeInstructions(fCompl, MODE, 1);

    fflush(fCompl);
    float* dEta_arr = safeMAlloc(MEASURES * sizeof(float));
    if (SCALE) {
        geomProgression(dEta_arr, (float)DETAMIN, (float)DETAMAX, MEASURES);
    } else {
        linearProgression(dEta_arr, (float)DETAMIN, (float)DETAMAX, MEASURES);
    }
    for (int i = 0; i < MEASURES; i++) {
        dEta = (double)dEta_arr[MEASURES - i - 1];
        int l = 0;

        char name[25];
        sprintf(name, "reinforcement_%0.5f", dEta);
        openFile(&fp1, name);
        writeInstructions(fp1, MODE, 0);
        buildInitialConfiguration();
        timeT = 0.0;
        evolveSystemTo(timeT + timeForWStat(dEta));
        while (l < LOOPS) {
            double t = timeT + correlationTime(dEta);
            evolveSystemTo(t);
            take2DMeasures(t);
            l++;
        }
        fprintf(fCompl, "%.5f", dEta);
        printGenerics(fCompl, 1);
        fclose(fp1);
    }
    fclose(fCompl);
    free(dEta_arr);
}

void onBorderWalkEvolution(void) {
    if (DIM != 1) {
        printf("Wrong mode for given dimension!");
        return;
    }
    float* deltaEtaArr = safeMAlloc((int)MEASURES * sizeof(float));
    if (SCALE) {
        geomProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, (int)MEASURES);
    } else {
        linearProgression(deltaEtaArr, (float)DETAMIN, (float)DETAMAX, (int)MEASURES);
    }
    openFile(&fCompl, "on_border_walk_MEAN");
    writeInstructions(fCompl, 1, 1);
    fprintf(fCompl, "# dEta eta0 eta1 eta2 W\n");
    for (int i = 0; i < (int)MEASURES; i++) {
        timeT = 0;
        dEta = deltaEtaArr[i];
        gen0 = 0.0, gen1 = 0.0, gen2 = 0.0, gen3 = 0.0;
        openFile(&fp1, "on_border_walk");
        writeInstructions(fp1, 1, 0);
        fprintf(fp1, "# eta0 eta1 eta2 W\n");
        buildInitialConfiguration();
        double nextTime = timeForWStat(dEta);
        evolveSystemTo(nextTime);
        int b1, b2, b1Last, b2Last;
        w1DBorders(&b1, &b2);
        b1Last = b1;
        b2Last = b2;
        int db1, db2;
        int l = 0;
        while (l < (int)LOOPS) {
            nextTime = timeT + 1. / (double)mobileCount;
            evolveSystemTo((double)nextTime);
            w1DBorders(&b1, &b2);
            db1 = b1 - b1Last; // compute the displacement based on last flip
            db2 = b2Last - b2;
            if (db1 > 0) {
                l++;
                gen0 += etaOfBiased(b1 + 1, b1), gen1 += etaOfBiased(b1 + 2, b1),
                    gen2 += etaOfBiased(b1 + 3, b1), gen3 += abs(b1 - b2);
                fprintf(fp1, "%.5f %.5f %.5f %d\n", etaOfBiased(b1 + 1, b1),
                        etaOfBiased(b1 + 2, b1), etaOfBiased(b1 + 3, b1), abs(b1 - b2));
                fflush(fp1);
            }
            if (db2 > 0) {
                l++;
                gen0 += etaOfBiased(b2 - 1, b2), gen1 += etaOfBiased(b2 - 2, b2),
                    gen2 += etaOfBiased(b2 - 3, b2), gen3 += abs(b1 - b2);
                fprintf(fp1, "%.5f %.5f %.5f %d\n", etaOfBiased(b2 - 1, b2),
                        etaOfBiased(b2 - 2, b2), etaOfBiased(b2 - 3, b2), abs(b1 - b2));
                fflush(fp1);
            }
            b1Last = b1;
            b2Last = b2;
        }
        fprintf(fCompl, "%.5f %.5f %.5f %.5f %.5f\n", dEta, gen0 / (double)LOOPS,
                gen1 / (double)LOOPS, gen2 / (double)LOOPS, gen3 / (double)LOOPS);
        fflush(fCompl);
        fclose(fp1);
    }
    free(deltaEtaArr);
    fclose(fCompl);
}

double etaOfBiased(int i0, int iBias) {
    if (s[i0] != s[iBias])
        return 0.0;
    return eta[i0];
}

double timeForWStat(double dEta) {
    return 1.0e2 / dEta;
}
double correlationTime(double dEta) {
    return 1.0e1 / dEta;
}

void evolveSystemTo(double nextTime) {
    while (timeT < nextTime) {
        if (mobileCount == 0) {
            mergeOldValues();
            break;
        }
        timeT += 1. / (double)mobileCount;
        if (timeT >= nextTime) {
            mergeOldValues();
        }
        singleInteraction();
    }
}

void inversionTime1D() {
    if (DIM != 1) {
        printf("Wrong mode!\n");
        return;
    }
    dEta = (double)DETA;
    char a[100] = "";
    sprintf(a, "step_time_dEtae%d", (int)log10(dEta));
    openFile(&fp1, a);
    writeInstructions(fp1, 3, 0);
    timeT = 0.0;
    double nTime, lastTime;
    buildInitialConfiguration();

    evolveSystemTo((double)timeForWStat(dEta)); //
    for (int i = 0; i < MEASURES; i++) {
        lastTime = timeT;
        int b01, b02, b1, b2, b1Last, b2Last;
        w1DBorders(&b01, &b02);
        b1 = b01;
        b1Last = b1;
        b2 = b02;
        b2Last = b2;
        int d1 = b1 - b1Last, d2 = b2 - b2Last;
        double invTime1 = 0, invTime2 = 0;
        int step1 = 0, step2 = 0;
        int f1 = 1, f2 = 1;
        while (f1 || f2) {
            nTime = timeT + 1. / (double)mobileCount;
            evolveSystemTo((double)nTime);
            w1DBorders(&b1, &b2);

            d1 = b1 - b1Last;
            d2 = b2 - b2Last;
            if (d1 < 0) {
                invTime1 = timeT - lastTime;
                step1 = b1 - b01 - d1;
                f1 = 0;
            }
            if (d2 < 0) {
                invTime2 = timeT - lastTime;
                step2 = b2 - b02 - d2;
                f2 = 0;
            }
            b1Last = b1;
            b2Last = b2;
        }
        fprintf(fp1, "%.2f %d %.2f %d %.2f %.2f\n", invTime1, step1, invTime2, step2,
                (float)(invTime1 + invTime2) / 2.0, (float)((step1 + step2) / 2.0));
        evolveSystemTo(timeT + (double)correlationTime(dEta));
    }
    fclose(fp1);
}

void w1DBorders(int* b1, int* b2) {
    int f1 = 1, f2 = f1, i1 = 0, i2 = lattice.dimension.size - 1;
    for (int i = 0; i < lattice.dimension.size; i++) {
        if (f1) {
            if (s[i] == s[i1] && z[i]) {
                *b1 = i;
            } else {
                f1 = 0;
            }
        }
        int j = i2 - i;
        if (f2) {
            if (s[j] == s[i2] && z[j]) {
                *b2 = j;
            } else {
                f2 = 0;
            }
        }
    }
}

void take1DMeasures(double nextTime) {
    switch (MODE) {
        case 0: {
            initiateGenerics(4);
            nZ = 0;
            nS1 = 0;
            nZ1 = 0;
            int lastZ0 = 0, lastZ1 = LSIZE - 1;
            for (int i = 0; i < N; i++) {
                nS1 += s[i];
                nZ1 += s[i] * z[i];
                nZ += z[i];
                if (z[i] && s[i] == s[0]) { //
                    lastZ0 = i;
                }
                int j = (int)LSIZE - i - 1;
                if (z[j] && s[j] == s[(int)LSIZE - 1]) {
                    lastZ1 = j;
                }
            }
            nClean = abs(lastZ1 - lastZ0) - 1;
            generics[0] += nS1;
            generics[1] += nZ1;
            generics[2] += nZ;
            generics[3] += nClean;
            fprintf(fp1, "%.2f %d %d %d %d\n", nextTime, nS1, nZ1, nZ, nClean);

            break;
        }
    }
}

void testtt(double t) {
    initiateGenerics(6);
    int* matrix_S0Z = safeMAlloc((int)N * sizeof(int));
    int* plot = safeMAlloc((int)N * sizeof(int));

    int* labSZ = safeMAlloc(N * sizeof(int)); // H-K label matrix of S;
    int* labS = safeMAlloc(N * sizeof(int));

    int* tempS = safeMAlloc((int)N * sizeof(int));
    int* tempZ = safeMAlloc((int)N * sizeof(int));

    int* onInnerLength = safeMAlloc(2 * N * sizeof(int));
    int* onOuterLength = safeMAlloc(2 * N * sizeof(int));

    double* maxHSPerX = safeMAlloc((int)(LSIZE - 2) * sizeof(double));
    double* maxHSZPerX = safeMAlloc((int)(LSIZE - 2) * sizeof(double));

    int inCount = 0, outCount = 0;
    for (int i = 0; i < N; i++) {
        plot[i] = s[i] + 2 * z[i];
        matrix_S0Z[i] = s[i] + 2 * z[i] == 2 ? 2 : -3;
        onInnerLength[i] = -1;
        onOuterLength[i] = -1;
        if (i < LSIZE - 2) {
            maxHSPerX[i] = NAN;
            maxHSZPerX[i] = NAN;
        }
    }

    labelCluster(labSZ, matrix_S0Z, &lattice);
    labelCluster(labS, s, &lattice);

    for (int i = 0; i < N; i++) {
        tempS[i] = labS[i] == labS[0];
        tempZ[i] = labSZ[i] == labSZ[0];
    }

    labelCluster(labS, tempS, &lattice);
    labelCluster(labSZ, tempZ, &lattice);

    int l_Inner = 0, l_Outer = 0;

    double hSM = 0, hSZM = 0;
    int hSCount = 0, hSZCount = 0;
    for (int i = 0; i < N; i++) {
        int x = i % LSIZE, y = (int)floor(i / (float)LSIZE);

        if (x == 0 || x == LSIZE - 1)
            continue;
        if (labS[i] == labS[0]) {
            for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
                int neig = lattice.neighbours->matrix[i][j];
                if (labS[neig] == labS[N - 1]) {
                    if (isValidInteraction(&lattice, i, neig)) {
                        plot[i] = -1;
                        l_Inner++;
                        if (j == 2 || j == 3) {
                            onInnerLength[inCount] = i;
                            inCount++;
                            plot[i] = -3;
                            hSM += y;
                            hSCount++;
                        }
                    }
                }
            }
        }

        if (labSZ[i] == labSZ[0]) {
            for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
                int neig = lattice.neighbours->matrix[i][j];
                if (labSZ[neig] == labSZ[N - 1]) {
                    if (isValidInteraction(&lattice, i, neig)) {
                        plot[i] = -2;
                        l_Outer++;
                        if (j == 2 || j == 3) {
                            onOuterLength[outCount] = i;
                            outCount++;
                            hSZM += y;
                            hSZCount++;
                        }
                    }
                }
            }
        }
    }

    hSM /= hSCount;
    hSZM /= hSZCount;

    double hS2M = 0, hSZ2M = 0;
    hSCount = 0, hSZCount = 0;
    for (int index = 0; index < inCount; index++) {
        int i = onInnerLength[index];
        int x = i % LSIZE, y = (int)floor(i / (float)LSIZE);
        bool f = maxHSPerX[x - 1] != maxHSPerX[x - 1];

        if (f || fabs(maxHSPerX[x - 1] - hSM) < fabs(y - hSM)) {
            maxHSPerX[x - 1] = y - hSM;
        }
    }

    for (int index = 0; index < outCount; index++) {
        int i = onOuterLength[index];
        int x = i % LSIZE, y = (int)floor(i / (float)LSIZE);
        bool f = maxHSZPerX[x - 1] != maxHSZPerX[x - 1];
        if (f || fabs(maxHSZPerX[x - 1]) < fabs(y - hSZM)) {
            maxHSZPerX[x - 1] = y - hSZM;
        }
    }

    for (int i = 0; i < LSIZE - 2; i++) {
        if (maxHSPerX[i] == maxHSPerX[i]) {
            hS2M += (double)maxHSPerX[i] * (double)maxHSPerX[i];
        }
        if (maxHSZPerX[i] == maxHSZPerX[i]) {
            hSZ2M += (double)maxHSZPerX[i] * (double)maxHSZPerX[i];
        }
    }

    hS2M /= (LSIZE - 2);
    hSZ2M /= (LSIZE - 2);

    if (EXPORT_MATRIX) {
        // exportMatrix(plot, N);
    }

    double w_Inner2 = hS2M;
    double w_Outer2 = hSZ2M;

    generics[0] = l_Inner;
    generics[1] = (l_Inner * l_Inner);
    generics[2] = w_Inner2;
    generics[3] = l_Outer;
    generics[4] = (l_Outer * l_Outer);
    generics[5] = w_Outer2;

    fprintf(fp1, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", t, (double)l_Inner, pow(l_Inner, 2),
            w_Inner2, (double)l_Outer, pow(l_Outer, 2), w_Outer2);
    fflush(fp1);
    free(matrix_S0Z);
    free(plot);
    free(labSZ);
    free(labS);
    free(tempS);
    free(tempZ);
    free(onOuterLength);
    free(onInnerLength);
    free(maxHSPerX);
    free(maxHSZPerX);
}

void take2DMeasures(double nextTime) {
    switch (MODE) {
        case 0:
            initiateGenerics(3);
            nZ = 0;
            nS1 = 0;
            nZ1 = 0;
            for (int i = 0; i < N; i++) {
                nS1 += s[i];
                nZ += z[i];
                nZ1 += s[i] * z[i];
            }

            generics[0] += nS1;
            generics[1] += nZ1;
            generics[2] += nZ;
            fprintf(fp1, "%.2f %d %d %d\n", nextTime, nS1, nZ1, nZ);

            break;
        case 1: {
            testtt(nextTime);
            break;
            initiateGenerics(6);
            int* matrix_SZ = safeMAlloc((int)N * sizeof(int));
            for (int i = 0; i < N; i++) {
                matrix_SZ[i] = s[i] + 2 * z[i];
            }

            int* labSZ = safeMAlloc(N * sizeof(int)); // H-K label matrix of S;
            int* labS = safeMAlloc(N * sizeof(int));
            labelCluster(labSZ, matrix_SZ, &lattice);
            labelCluster(labS, s, &lattice);

            int l_Inner = 0, inCount = 0;
            double innerY = 0, innerY2 = 0;

            int l_Outer = 0, outCount = 0;
            double outerY = 0, outerY2 = 0;

            for (int i = 0; i < N; i++) {
                int x = i % LSIZE, y = (int)floor(i / (float)LSIZE);
                if (x == 0 || x == LSIZE - 1)
                    continue;
                bool in = false, out = false;
                if (labS[i] == labS[0]) {
                    for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
                        if (labS[lattice.neighbours->matrix[i][j]] != labS[0]) {
                            if (isValidInteraction(&lattice, i, lattice.neighbours->matrix[i][j])) {
                                l_Inner++;
                                in = true;
                            }
                        }
                    }
                }

                if (labSZ[i] == labSZ[0]) {
                    for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
                        if (labSZ[lattice.neighbours->matrix[i][j]] != labSZ[0]) {
                            if (isValidInteraction(&lattice, i, lattice.neighbours->matrix[i][j])) {
                                l_Outer++;
                                out = true;
                            }
                        }
                    }
                }
                if (in) {
                    inCount++;
                    innerY += y;
                    innerY2 += y * y;
                }
                if (out) {
                    outCount++;
                    outerY += y;
                    outerY2 += y * y;
                }
            }
            if (inCount > 0) {
                innerY /= (float)inCount;
                innerY2 /= (float)inCount;
            }
            if (outCount > 0) {
                outerY /= (float)outCount;
                outerY2 /= (float)outCount;
            }
            double w_Inner2 = innerY2 - (innerY * innerY);
            double w_Outer2 = outerY2 - (outerY * outerY);

            generics[0] = l_Inner;
            generics[1] = (l_Inner * l_Inner);
            generics[2] = w_Inner2;
            generics[3] = l_Outer;
            generics[4] = (l_Outer * l_Outer);
            generics[5] = w_Outer2;

            fprintf(fp1, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", nextTime, (double)l_Inner,
                    pow(l_Inner, 2), w_Inner2, (double)l_Outer, pow(l_Outer, 2), w_Outer2);
            /*
            double w_Inner = 0, w_Outer1 = 0;
            int l_Inner = 0, l_Outer1 = 0;
            l_Inner = mc_biasedwalk((LSIZE / 2 - 1) * LSIZE, (N / 2) - 1, labS, &w_Inner);
            l_Outer1 = mc_biasedwalk((LSIZE / 2 - 1) * LSIZE, (N / 2) - 1, labSZ, &w_Outer1);

            generics[0] = l_Inner;
            generics[1] = pow(l_Inner, 2);
            generics[2] = w_Inner;
            generics[3] = l_Outer1;
            generics[4] = pow(l_Outer1, 2);
            generics[5] = w_Outer1;

            fprintf(fp1, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", nextTime, (double)l_Inner,
                    pow(l_Inner, 2), w_Inner, (double)l_Outer1, pow(l_Outer1, 2), w_Outer1);
            */
            fflush(fp1);
            free(matrix_SZ);
            free(labSZ);
            free(labS);
            break;
        }
        case 2: {
            initiateGenerics(6);
            double MZ0 = 0, MZ1 = 0, MS0 = 0, M2Z0 = 0, M2Z1 = 0, M2S0 = 0;
            for (int x = 1; x < LSIZE - 1; x++) {
                int lastZ0 = -1, lastZ1 = -1, lastS0 = -1;
                for (int y = 0; y < LSIZE; y++) {
                    int i = x + LSIZE * y;
                    if (lastZ1 == -1) {
                        if (s[i] == 1 && z[i] == 1) {
                            lastZ1 = y;
                        }
                    }
                    if (s[i] == 0) {
                        lastS0 = y;
                        if (z[i] == 1) {
                            lastZ0 = y;
                        }
                    }
                }

                MZ0 += lastZ0;
                MZ1 += lastZ1;
                MS0 += lastS0;
                M2Z0 += pow(lastZ0, 2);
                M2Z1 += pow(lastZ1, 2);
                M2S0 += pow(lastS0, 2);
            }
            MZ0 /= (LSIZE - 2);
            MZ1 /= (LSIZE - 2);
            MS0 /= (LSIZE - 2);
            M2Z0 /= (LSIZE - 2);
            M2Z1 /= (LSIZE - 2);
            M2S0 /= (LSIZE - 2);

            fprintf(fp1, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", nextTime, MS0, MZ0, MZ1, M2S0,
                    M2Z0, M2Z1);
            fflush(fp1);

            generics[0] += MS0;
            generics[1] += MZ0;
            generics[2] += MZ1;
            generics[3] += M2S0;
            generics[4] += M2Z0;
            generics[5] += M2Z1;
        }
    }
}

void initiateGenerics(int length) {
    if (length <= 0 || generics != NULL) {
        return;
    }
    genericsLength = length;
    generics = safeMAlloc(genericsLength * sizeof(double));
    for (int i = 0; i < genericsLength; i++) {
        generics[i] = 0;
    }
}

void printGenerics(FILE* f, int shouldFree) {
    if (fCompl == NULL || generics == NULL || genericsLength <= 0) {
        printf("WARNING: trying to write generic measures array with invalid structure");
        return;
    }
    for (int i = 0; i < genericsLength; i++) {
        fprintf(f, " %.5f", generics[i] / LOOPS);
    }
    fprintf(f, "\n");
    fflush(f);
    if (shouldFree) {
        free(generics);
        generics = NULL;
    }
}

void takeMeasures(double nextTime) {

    if (DIM == 1) {
        take1DMeasures(nextTime);
    } else {
        take2DMeasures(nextTime);
    }
    fflush(fp1);
    if (EXPORT_MATRIX)
        exportBaseSpinMatrix();
}

void clear(void) {
    free(s);
    free(post_s);
    free(z);
    free(post_z);
    free(eta);
    free(post_eta);
    free(mobile);
    free(locMobile);
    destroyLattice(&lattice);
}

//==================Virtual Environment==================//
//      It acts wit post_s and post_z, virtual versions
//      of s and z matrixes. When the time of the simu-
//      lation gets closer to the time for the measure
//      the virtual values from post_s and post_z are
//      merged for s and z. Never use post_s and post_z
//      to take measures, values from this matrixes do not
//      correspond for the values of timeT.
void singleInteraction(void) {
    int changedState = 0;

    int pPos = mobile[randomInt(mobileCount)];
    int vPos = lattice.neighbours->matrix[pPos][randomInt(lattice.neighbours->neighboursCount)];

    if (!isValidInteraction(&lattice, pPos, vPos)) {
        // printf("%d %d\n", pPos, vPos);
        return;
    }
    if (FIXED_BORDERS) {
        int on1DBorder = (pPos == 0 || pPos == LSIZE - 1);
        if (DIM == 1 && on1DBorder) {
            if (C == 2) {
                if (pPos == 0)
                    return;
            } else {
                return;
            }

        } else if (DIM == 2) {
            int x = pPos % LSIZE, y = pPos / LSIZE;
            int on2DBorder = y == 0 || y == HSIZE - 1;
            on2DBorder |= (B == 3 && (x == 0 || x == LSIZE - 1));
            if (on2DBorder) {
                // printf("%d\n", pPos);
                return;
            }
        }
    }
    int reinforcementInteraction = post_s[pPos] == post_s[vPos];

    if (reinforcementInteraction) {
        post_eta[pPos] += (double)dEta;
        if (post_eta[pPos] >= 1.0) {
            post_eta[pPos] = 1.0;
            post_z[pPos] = 1;
            changedState = 1;
        }
    } else {
        if (post_z[pPos] == 1) {
            post_z[pPos] = 0;
            post_eta[pPos] = 0.0;
        } else {
            post_s[pPos] = post_s[vPos];
            post_eta[pPos] = 0.0;
            changedState = 1;
        }
    }

    if (changedState) {
        updateMobileMatrix(pPos);
        for (int i = 0; i < lattice.neighbours->neighboursCount; i++) {
            updateMobileMatrix(lattice.neighbours->matrix[pPos][i]);
        }
    }
    return;
}
void updateMobileMatrix(int i) {

    if (isMobile(i)) {

        if (locMobile[i] >= mobileCount) { // is mobile but its not on the mobile list

            // change pos of i with the last element of the mobile

            locMatrixUpdate(i, mobile[mobileCount], mobile, locMobile);
            mobileCount++;
        }
    } else {
        if (locMobile[i] < mobileCount) { // its not mobile but is on the llist
            // change pos of i with the last element of the mobile
            locMatrixUpdate(i, mobile[mobileCount - 1], mobile, locMobile);

            mobileCount--;
        };
    }
}
int isMobile(int i) {
    bool f = dEta == 0.0;
    if (post_z[i] == 1 || f) {
        for (int j = 0; j < lattice.neighbours->neighboursCount; j++) {
            int neigIndex = lattice.neighbours->matrix[i][j];
            if (post_s[neigIndex] != post_s[i]) {
                return 1;
            }
        }
        return 0;
    }
    return 1;
}
void locMatrixUpdate(int s1, int s2, int* mat, int* locMat) {
    int tempS;

    mat[locMat[s1]] = s2;
    mat[locMat[s2]] = s1;

    tempS = locMat[s1];
    locMat[s1] = locMat[s2];
    locMat[s2] = tempS;
    return;
}

void openFile(FILE** f, char* prefix) {
    if (seed == -1) {
        seed = setupRandom(SEED);
    }
    char name[120];
    sprintf(name, "data_pvm_%s", prefix);
    *f = safeSeedOpen(name, ".dat", &seed, SEED != 0);
    fflush(*f);
}

void writeInstructions(FILE* f, int mode, int isCompl) {
    writeHeader(f, "Persistent Voter Model: pvm.c", "L. Possamai");
    fprintf(f, "# ║ Seed: %d\n", seed);
    fprintf(f, "# ║ Dimension: %d\n", DIM);

    if (SCALE) {
        fprintf(f, "# ║ Log Measures: %d\n", (int)MEASURES);
    } else {
        fprintf(f, "# ║ Linear Measures: %d\n", (int)MEASURES);
    }
    switch (EVOLUTION) {
        case 0:
            fprintf(f, "# ║ dEta = %.5f\n", DETA);
            fprintf(f, "# ║ L = %d\n", (int)LSIZE);
            fprintf(f, "# ║ H = %d\n", (int)HSIZE);
            fprintf(f, "# ║ t: [%d, %d]\n", (int)TMIN, (int)TMAX);
            break;
        case 1:
            fprintf(f, "# ║ L = %d\n", (int)LSIZE);
            fprintf(f, "# ║ H = %d\n", (int)HSIZE);
            if (isCompl) {
                fprintf(f, "# ║ loops: %d\n", (int)LOOPS);
                fprintf(f, "# ║ dEta: [%.5f, %.5f]\n", (double)DETAMIN, (double)DETAMAX);
            } else {
                fprintf(f, "# ║ dEta: %.5f\n", (double)dEta);
            }
            break;
        case 2:
            fprintf(f, "# ║ dEta = %.5f\n", DETA);
            fprintf(f, "# ║ H = %d\n", (int)HSIZE);
            if (isCompl) {
                fprintf(f, "# ║ loops: %d\n", (int)LOOPS);
                fprintf(f, "# ║ L: [%d, %d]\n", (int)LMIN, (int)LMAX);
            } else {
                fprintf(f, "# ║ L: %d\n", (int)lattice.dimension.proportions[0]);
            }
            break;
        case 8:
        case 9:
            fprintf(f, "# L = %d\n", (int)LSIZE);
            if (isCompl) {
                fprintf(f, "# loops: %d\n", (int)LOOPS);
                fprintf(f, "# dEta: [%.5f, %.5f]\n", DETAMIN, (double)DETAMAX);
            } else {
                fprintf(f, "# dEta: %.5f\n", (double)dEta);
            }
            break;
    }
    if (DIM == 1) {
        write1DInstructions(f, MODE, isCompl);
    } else if (DIM == 2) {
        write2DInstructions(f, MODE, isCompl);
    }
    fflush(f);
}

void write2DInstructions(FILE* f, int mode, int isCompl) {
    switch (mode) {
        case 0:
            fprintf(f, "# ╚ t nS0 nZ0 nZ\n");
            break;
        case 1:
            fprintf(f, "# ╚ t | lInner | lInner2 | wInner | lOuter | lOuter2 | wOuter \n");
            break;
        case 2:
            if (isCompl) {
                fprintf(f, "# ╚ dEta <yS0> <yZ0> <yZ1> <y2S0> <y2Z0> <y2Z1> \n");
            } else {
                fprintf(f, "# ╚ t <yS0> <yZ0> <yZ1> <y2S0> <y2Z0> <y2Z1> \n");
            }
            break;
    }
}
void write1DInstructions(FILE* f, int mode, int isCompl) {

    switch (mode) {
        case 0:
            if (isCompl) {
                fprintf(f, "# ╚ dEta N_S0 N_Z0 N_Z N_free\n");
            } else {
                fprintf(f, "# ╚ t N_S0 N_Z0 N_Z N_free\n");
            }
            break;
        case 1:
            break;
        case 3:
            fprintf(f, "# ╚ t step\n");
        case 8:
            if (isCompl) {
                fprintf(f, "# ╚ dEta <vel>\n");
            } else {
                fprintf(f, "# ╚ t vel\n");
            }
            break;
        case 9:
            if (isCompl) {
                fprintf(f, "#  dEta <eta0>\n");
            } else {
                fprintf(f, "#  t eta0\n");
            }
    }
}

/*************************************************************************
 * *                     Biased walks along the external hull               *
 * *                          Last modified: 17/06/2024                     *
 * * Returns the area enclosed by the hull, despite the presence of smaller *
 * * internal domains. We depart from the cluster labelling site (smaller   *
 * * index), and the initial contribution to the area is LSIZE+1. Then we   *
 * * attempt to walk along the left, front, right or backward directions.   *
 * * The incoming (backward) direction is only accepted if it's not         *
 * * possible to go on the other directions. We choose the height of the    *
 * * initial size as 10*LSIZE and update the height along the walk.         *
 * * Percolating clusters are not taken into account since the walls are    *
 * * disconnected in that case. The area is also updated along the walk     *
 * * using  the following table (the column shows the precedent step):      *
 * *                                                                        *
 * *              up     right    down   left                               *
 * *  up (0)       0      h+1      h+1     0                                *
 * *  right (3)    0      h+1      -h      1                                *
 * *  down (2)    -h      0        0      -h                                *
 * *  left (1)    -h      1        0      h+1                               *
 * *                                                                        *
 * *                                                                        *
 * *************************************************************************/

int mc_biasedwalk(int qual, int end, int* lab, double* w2) {
    int i = 0, dir = 0, ok;

    if (lab[qual] != lab[end]) {
        printf("Unable to calculate length of interface, start and end points are not connected");
        return -1;
    }

    /* First check whether the starting point has one or two branches going out. If there
     *    are two branches connected by the starting point, the configuration shoulbe be:  11
     *                                                                                        10
     *                                                                                           The
     * walk will first go through the horizontal branch and, in order to go to the down branch
     * it should pass over the initial site, but we should not stop the walk there. To do this,
     * we set endpoint=0. When the walk returns from the horizontal branch, we set it to 1, so
     * the walk can stop the next time it visits the starting point. Notice however that the
     * horizontal branch may close the loop and join the down one without passing through the
     * starting point. For example:  111 101 111                              */
    int iR = lattice.neighbours->matrix[qual][0], iL = lattice.neighbours->matrix[qual][1],
        iD = lattice.neighbours->matrix[qual][2], iU = lattice.neighbours->matrix[qual][3];

    /* from the starting point, choose the direction to move, there are just 2 possibilities
     * (from the way we choose it): */
    if (lab[iR] == lab[qual]) {
        dir = 3;
        i = iR;
    } else if (lab[iD] == lab[qual]) {
        dir = 2;
        i = iD;
    }

    int* arr = safeMAlloc(N * sizeof(int));

    for (int i = 0; i < N; i++) {
        arr[i] = s[i] + 2 * z[i];
    }
    /* start the walk around the cluster, clockwise: */

    int length = 0, count = 0;
    double yM = 0, y2M = 0;
    while (i != end)
    /* while the walk doesn't return to the starting point, or if it returns, it can continue to
       the other branch */
    {
        if (i == end)
            break;

        dir = (dir + 1) % 4; /* from the incoming direction, try left first */

        iR = lattice.neighbours->matrix[i][0], iL = lattice.neighbours->matrix[i][1],
        iD = lattice.neighbours->matrix[i][2], iU = lattice.neighbours->matrix[i][3];

        int l = (lab[iR] != lab[qual]) + (lab[iU] != lab[qual]) + (lab[iL] != lab[qual]) +
                (lab[iD] != lab[qual]);

        int x = i % LSIZE;
        int yi = (int)floor(i / (float)LSIZE);
        if (x != 0 && x != LSIZE - 1) {
            length += l;
            yM += yi;
            y2M += pow(yi, 2);
            count++;
        }
        if (l != 0) {
            arr[i] = -2;
        }

        ok = 0;
        while (!ok) /* from the incoming direction: try right, in front, left and backwards */
        {
            switch (dir) {
                case 0:
                    if (lab[iU] == lab[i]) {
                        ok = 1;
                        i = iU;
                    }
                    break;
                case 1:
                    if (lab[iL] == lab[i]) {
                        ok = 1;
                        i = iL;
                    }
                    break;
                case 2:
                    if (lab[iD] == lab[i]) {
                        ok = 1;
                        i = iD;
                    }
                    break;
                case 3:
                    if (lab[iR] == lab[i]) {
                        ok = 1;
                        i = iR;
                    }
                    break;
            }
            if (ok == 0)
                dir = (dir + 3) % 4;
        }
    }
    yM /= count;
    y2M /= count;
    *w2 = y2M - pow(yM, 2);
    if (EXPORT_MATRIX)
        exportMatrix(arr, N);
    free(arr);
    return length;
}
