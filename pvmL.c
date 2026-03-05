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

//==================Lattice Config.==================//

#define SEED 0 // 123456789: random

//==================Lattice Config.==================//

#define LSIZE 150
#define DIM 2

//==================Generic.Evolution================//

#define DETA 1e-3

//==================Temp. Evolution==================//

#define TMIN 1e0 // cant be 0 for scale log
#define TMAX 1e5

//==================Length Evolution=================//

#define LMAX 256
#define LMIN 8

//==================DEta. Evolution==================//

#define DETAMIN 1e-4
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
#define MEASURES 20

// Number of loops(only valid for non temporal evolutions)
#define LOOPS 1e3

// 0: Linear
// 1: LOG
#define SCALE 1

// 0: Number of zealots
// 1: Interface Length
// 2: Instantaneous AR Width
#define MODE 1

/*
 * 0: temporalEvolution
 * 1: reinforcementEvolution
 * 2: lengthEvolution
 * 3: inversionTime1D
 * 4: onBorderWalkEvolution
 */
#define EVOLUTION 0

//=================Boundry Condition=================//
#define FIXED_BORDERS 1

/*
 * "b" = "boundry condition"
 * 0 -> periodic
 * 1 -> dobrushin vertical(periodic horizontally)
 * 2 -> dobrushin horizontal(periodic vertically)
 * 3 -> square
 */
#define B 3
#define EXPORT_MATRIX 0

// ==================================================//

int N;
SquareLattice lattice;
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
    N = pow(LSIZE, DIM);
    initiateSquareLattice(&lattice, DIM, LSIZE, B);

    if (EXPORT_MATRIX)
        fMatrix = fopen("matrix.dat", "w");
    setup();
    initiate();
    clear();
    if (EXPORT_MATRIX)
        fclose(fMatrix);
    return 0;
}

void initiate(void) {

    switch (EVOLUTION) {
        case 0: // temporalEvolution
            temporalEvolution();
            fclose(fp1);
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

void buildInitialConfiguration(void) {
    switch (C) {
        case 0:
            for (int i = 0; i < N; i++) {
                post_s[i] = randomInt(2);
                post_z[i] = 1;
                post_eta[i] = 1.0;
            }
            break;
        case 1:
            for (int i = 0; i < N; i++) {
                post_s[i] = i < (N / 2) ? 0 : 1;
                post_z[i] = 1;
                post_eta[i] = 1.0;
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
                post_s[i] = distance <= (float)lattice.L / 4.0f ? 0 : 1;
                post_z[i] = 1;
                post_eta[i] = 1.0;
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

void exportMatrix(int* s, int sizeN) {
    for (int i = 0; i < sizeN; i++) {
        fprintf(fMatrix, "%d ", s[i]);
    }
    fprintf(fMatrix, "\n");
}

void exportBaseSpinMatrix() {
    int* a = safeMAlloc(N * sizeof(int));
    for (int i = 0; i < N; i++) {
        a[i] = s[i] + 2 * z[i];
    }
    int* b = safeMAlloc(N * sizeof(int));
    labelClusterEff(b, s, &lattice);
    if (0) {
        exportMatrix(a, N);
    } else {
        exportMatrix(b, N);
    }
    free(a);
    free(b);
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
    buildInitialConfiguration();
    dEta = (double)DETA;

    if (EXPORT_MATRIX)
        exportBaseSpinMatrix();
    for (int i = 0; i < MEASURES; i++) {
        nextTime = (double)time_arr[i];
        evolveSystemTo(nextTime);
        takeMeasures(nextTime);
    }
    free(time_arr);
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
        initiateSquareLattice(&lattice, DIM, length, B);
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
    int f1 = 1, f2 = f1, i1 = 0, i2 = lattice.size - 1;
    for (int i = 0; i < lattice.size; i++) {
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
            initiateGenerics(2);
            int* matrix_SZ = safeMAlloc((int)N * sizeof(int));
            for (int i = 0; i < N; i++) {
                matrix_SZ[i] = s[i] + 2 * z[i];
            }

            int* labSZ = safeMAlloc(N * sizeof(int)); // H-K label matrix of S;
            int* labS = safeMAlloc(N * sizeof(int));
            labelClusterEff(labSZ, matrix_SZ, &lattice);
            labelClusterEff(labS, s, &lattice);

            int lOuter = 0, lInner = 0;

            for (int i = 0; i < N; i++) {

                if (labSZ[i] == labSZ[0] ||
                    labSZ[i] == labSZ[N - 1]) { // is from one of the bigger clusters
                    for (int j = 0; j < (int)lattice.neighbours->neighboursCount; j++) {
                        int iNeig = lattice.neighbours->matrix[i][j];
                        if (isValidInteraction(&lattice, i, iNeig) && labSZ[i] != labSZ[iNeig])
                            lOuter++;
                    }
                }
                if (labS[i] == labS[0]) { //
                    for (int j = 0; j < (int)lattice.neighbours->neighboursCount; j++) {
                        int iNeig = lattice.neighbours->matrix[i][j];
                        if (isValidInteraction(&lattice, i, iNeig) && s[i] != s[iNeig])
                            lInner++;
                    }
                }
            }

            generics[0] += (double)lOuter / 2;
            generics[1] += (double)lInner;
            fprintf(fp1, "%.2f %.2f %.2f\n", nextTime, (double)lOuter / 2.0, (double)lInner);
            fflush(fp1);
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
            int on2DBorder = y == 0 || y == LSIZE - 1;
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
    if (post_z[i] == 1) {
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
            fprintf(f, "# ║ t: [%d, %d]\n", (int)TMIN, (int)TMAX);
            break;
        case 1:
            fprintf(f, "# ║ L = %d\n", (int)LSIZE);
            if (isCompl) {
                fprintf(f, "# ║ loops: %d\n", (int)LOOPS);
                fprintf(f, "# ║ dEta: [%.5f, %.5f]\n", (double)DETAMIN, (double)DETAMAX);
            } else {
                fprintf(f, "# ║ dEta: %.5f\n", (double)dEta);
            }
            break;
        case 2:
            fprintf(f, "# ║ dEta = %.5f\n", DETA);
            if (isCompl) {
                fprintf(f, "# ║ loops: %d\n", (int)LOOPS);
                fprintf(f, "# ║ L: [%d, %d]\n", (int)LMIN, (int)LMAX);
            } else {
                fprintf(f, "# ║ L: %d\n", (int)lattice.L);
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
            fprintf(f, "# ╚ - | lOuter | lInner \n");
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
