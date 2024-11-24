#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mc2023.h>
#include <lat2eps.h>
#include <monte_carlo.h>


void teste(void);
void post(void);
void pre(void);
void consensus(void);
void evolucao(void);
void criarPropriedades(void);
void populacionaVizinhos(void);
void atribuirPropriedades(void);
void criarVizinhos(void);
int bloqueado(int);
void interacoes(void);
unsigned long generateFile(void);
void hoshen_kopelman(int*, int*, int*);
void verificarLiberdadeVizinhos(int, int);
int percolates2d(int, int*);
void createPerfectLogTable(int*);
void comecar(int);
void iniciar_coleta(int);
void verificarCruzamentoInterface(void);
void save_configuration(int, char*);
void time_evolution(void);
void write(double, int, int, int);
void updateOldValues(void);
void interface(void);
int getInterfacialSizeOfCluster(int*, int*);
void getBiggestCluster(int*, int*, int);
void unique_cluster(void);


/*******************************************************************/
/**                             CONSTANTES                        **/
/*******************************************************************/
#define LOOPS                   500
#define TEMPO_MAX               1e5
#define PASSO_TEMPO             500
#define DELTAETA_INICIAL        1./1000
#define LSIZE                   128
#define NSIZE                   LSIZE * LSIZE


#define EVOLUCAO                2           // 0: NORMAL
                                            // 1: ATINGIR CONSENSO
                                            // 2: Interface Lenght

#define CONFIGURACAO_I          2           // 1 : CIRCULAR
                                            // 2 : LINEAR**
                                            // 3 : ALEATORIO
/********** Par�metros para EVOLUCAO = 0(normal) *******************/
#define ESCALA                  1           // 0: LINEAR
                                            // 1: LOG
/********** Par�metros para ESCALA = 1(log) ************************/
#define MEDIDAS                 140//70

/********** Par�metros para EVOLUCAO = 1(consensus) ****************/
#define DELTAETA_FINAL          1.0
#define DELTAETA_MEDIDAS        41

/******* Par�metros para CONFIGURACAO_I = CIRCULAR *****************/
#define RAIO                    10          //NAO IMPLEMENTADO AINDA

#define FORCE_SEED              1
/******* Par�metros para CONDICAO DE CONTORNO LINEAR FIXA **********/
#define FIXED_EXTREMES          1           //CILINDRIC
#define FIXED_BORDERS           0           //SQUARED
#define L_INICIAL               16
#define L_FINAL                 256
/*******************************************************************/
#define SAVE_CONFIG             1
#define CONFIG_T0               0
#define CONFIG_DELTA_T          100000
#define CONFIG_TF               111000//CONFIG_T0 + 10 * CONFIG_DELTA_T



int *dir, *esq, *cima, *baixo, **vizinhos, *zelotes, *polaridade/*0, 1 */,
    *moveis, *quaisMoveis, *tempos, mag, ncl, zar, interfCol;
int nMoveis0, mag0, *polaridade0, *zelotes0,  *moveis0, *quaisMoveis0,
    l0_int = 0, l1_int = 0, L = LSIZE, N = NSIZE;
float *confianca0;
int *v, nMoveis;
double tempo,deltaEta;
float *confianca, *o_confianca;
unsigned int seed;
FILE *fp1;

int main(void){
    int lps = 200;//getIntInputFromClient("loops = ", 1, 1e4);
    if(FORCE_SEED){
        lps = 1;
    }
    iniciar_coleta(lps);

}

void iniciar_coleta(int i){
    int j, t = 0;
    printf("  Starting %d loops for data collection. \n", i);
    for(j = 0; j < i; j++){
        t = 0 - clock();
        comecar(j);
        t += clock();
        printf("  Collection concluded, duration time: %.2f seconds \n",(double) t/CLOCKS_PER_SEC);
    }
    return;
}


void comecar(int i){
    printf("\n    Starting the %dth collection\n", i + 1);
    pre();
    teste();
    if(EVOLUCAO == 0) {
        time_evolution();
    } else if(EVOLUCAO == 1){
        consensus();
    } else if(EVOLUCAO == 2){
        interface();
    }
    post();
}

void teste(void){

}



void post(void){
    free(dir);
    free(esq);
    free(cima);
    free(baixo);
    free(vizinhos);
    free(v);
    free(tempos);
    free(zelotes);
    free(confianca);
    free(polaridade);
    free(moveis);
    free(quaisMoveis);
    free(zelotes0);
    free(confianca0);
    free(polaridade0);
    free(moveis0);
    free(quaisMoveis0);

    fclose(fp1);
}



void pre(){
    seed=0;
    seed = generateFile();
    start_randomic(seed);

    criarVizinhos();
    criarPropriedades();


    create_time_table(tempos, TEMPO_MAX, (int) MEDIDAS - 1, ESCALA);

    return;
}
void criarPropriedades(void){
    zelotes = jmalloc(N * sizeof(int));
    zelotes0 = jmalloc(N * sizeof(int));
    confianca = jmalloc(N * sizeof(float));
    confianca0 = jmalloc(N * sizeof(float));
    polaridade = jmalloc(N * sizeof(int));
    polaridade0 = jmalloc(N * sizeof(int));
    moveis = jmalloc(N * sizeof(int));
    moveis0 = jmalloc(N * sizeof(int));
    quaisMoveis = jmalloc(N * sizeof(int));
    quaisMoveis0 = jmalloc(N * sizeof(int));

    tempos = jmalloc(sizeof(int) * MEDIDAS);

}

void updateOldValues(void){
    nMoveis0 = nMoveis;
    mag0 = mag;
    zelotes0 = zelotes;
    confianca0 = confianca;
    polaridade0 = polaridade;
    moveis0 = moveis;
    quaisMoveis0 = quaisMoveis;
    return;
}


void consensus(void){
    fprintf(fp1,"#  Deta   t  \n");
    deltaEta = DELTAETA_INICIAL;
    while(deltaEta <= DELTAETA_FINAL){
        atribuirPropriedades();
        tempo = 0.0;
        while(tempo < TEMPO_MAX){
            tempo += 1./nMoveis;
            interacoes();
            if((mag == 0) || (mag == N)){
                fprintf(fp1,"%.7f %.6f\n",deltaEta,tempo);
                fflush(fp1);
                break;
            }
            if(tempo > TEMPO_MAX){
                fprintf(fp1,"# %.7f %d MAX\n",deltaEta,mag);
                  fflush(fp1);
                  break;
            }
        }
        deltaEta *= pow(DELTAETA_FINAL / DELTAETA_INICIAL, 1. / (DELTAETA_MEDIDAS - 1.));
    }
}

void interface(void){
  int measures = 20;
  float *l_arr = jmalloc(20 * sizeof(float));
  geomProgression(l_arr, L_INICIAL, L_FINAL, measures);

  char name[50];
    L = L_INICIAL;
    int threshold = 1e3, l0, l1, lar;
    deltaEta = 1e-2;
    threshold = 1e3;
    fprintf(fp1,"#  deltaEta: %.5f\n", deltaEta);
    fprintf(fp1,"#  threshold: %d\n", threshold);
    fprintf(fp1,"#  L    <L0>   <L1>   <L_Ar>\n");
    for(int i = 0; i < measures; i++){
      L = (int) l_arr[i];
      N = L * L;
        sprintf(name, "%d", L);
        criarVizinhos();
        criarPropriedades();
        atribuirPropriedades();
        updateOldValues();
        tempo = 0.0;
        l0 = 0; l1 = 0; lar = 0;
        while(tempo < TEMPO_MAX){
          tempo += 1./nMoveis;
            interacoes();
            updateOldValues();
            if(tempo >= threshold){
              unique_cluster();
              l0 = l0_int;
              l1 = l1_int;
              fprintf(fp1,"%d %d %d %d\n", L, l0 - L, l1 - L, lar);
              fflush(fp1);

              if(SAVE_CONFIG) save_configuration(tempo, name);
              break;
            }
            if(tempo > TEMPO_MAX){
                break;
            }
        }
    }
}



void evolucao(void){
    int i, z = 0, z1 = 1, per, ti = 0, tiMax = (int) MEDIDAS - 1, oldTempo = CONFIG_T0;
    double proxTempo;
    deltaEta = DELTAETA_INICIAL;
    fprintf(fp1,"# Deta = %.5f\n",deltaEta);
    fprintf(fp1,"#  t  m  z1  z  ncl z_ar per mob int_cross\n");
    atribuirPropriedades();
    tempo = 0.0;
    while(ti <= tiMax){
        proxTempo = tempos[ti];
        while(tempo <= proxTempo && nMoveis > 0){
            interacoes();
            tempo += 1./nMoveis;
            if(nMoveis == 0) tempo = tempos[tiMax];

            if(SAVE_CONFIG == 1){
                if (tempo >= CONFIG_T0 && tempo <= CONFIG_TF) {
                    if(tempo > oldTempo){
                        printf("flash");
                        save_configuration(tempo, "a");
                    }
                }
            }
        }


        if(tempo >= proxTempo){
                z = 0;
                z1 = 0;
                for(i = 0; i < N; i++){
                    if(polaridade[i] == 1){
                        z1 += zelotes[i];
                    }
                    z += zelotes[i];
                }
                hoshen_kopelman(NULL,&ncl,&per);
                verificarCruzamentoInterface();

            fprintf(fp1,"%.5f  %d %d %d %d %d %d %d %d\n", proxTempo, mag, z1, z, ncl, zar, per, nMoveis, interfCol);
            fflush(fp1);
            do{
                ti += 1;
            } while(tempos[ti] <= proxTempo);


        }
    }
    return;
}

void write(double proxTempo, int z, int z1, int per){
    fprintf(fp1,"%.5f  %d %d %d %d %d %d %d %d\n", proxTempo, mag0, z1, z, ncl, zar, per, nMoveis0, interfCol);
    if(SAVE_CONFIG == 1){
        if (proxTempo >= CONFIG_T0 && proxTempo <= CONFIG_TF) {
            save_configuration(proxTempo, "evolution");
        }

    }
    fflush(fp1);
}

void time_evolution(void){
    double proxTempo = 0;
    int tiMax = MEDIDAS - 1, ti = 1, z=0, z1=1, per;
    deltaEta = DELTAETA_INICIAL;
    fprintf(fp1,"# Deta = %.5f\n",deltaEta);
    fprintf(fp1,"#  t  m  z1  z  ncl z_ar per mob int_cross\n");
    double nT = 0;
    atribuirPropriedades();
    updateOldValues();
    tempo = 0.0;
    proxTempo = tempos[ti];

    while(nMoveis > 0 && tempo <= tempos[tiMax]){
        interacoes();
        if(nMoveis == 0) {
            nT = tempos[tiMax];
        } else {
            nT = tempo + 1./nMoveis;
        }
        if(nT > proxTempo){
            z = 0;
            z1 = 0;
            for(int i = 0; i < N; i++){
                if(polaridade0[i] == 1){
                    z1 += zelotes0[i];
                }
                z += zelotes0[i];
            }
            hoshen_kopelman(NULL,&ncl,&per);
            verificarCruzamentoInterface();
            while(nT > proxTempo){
                write(proxTempo, z, z1, per);//OLDONES
                ti += 1;
                //printf("%.2f %.2f", proxTempo, nT);
                if(ti > tiMax) break;
                proxTempo = tempos[ti];
            }
        }
        tempo = nT;
        updateOldValues();
    }
    return;
}

void interacoes(void){
    int pPos,pIndice, vIndice, vPos, virouZ = 0, virouN = 0, pZelote, vZelote, pIgual;

    pIndice = (int) (FRANDOM * nMoveis);
    pPos = moveis[pIndice];
    vIndice = (int) (FRANDOM * 4);
    vPos = vizinhos[pPos][vIndice];

    if(CONFIGURACAO_I == 2){
        if((FIXED_EXTREMES == 1 || FIXED_BORDERS == 1) && (pPos < L || pPos >= N - L)) return;
        if(FIXED_BORDERS == 1 && (((pPos % L) == L - 1) || (pPos % L == 0))) return;
        if((vPos < L && vIndice == 3)|| (pPos < L && vIndice == 2)) return; //prevents top/bottom interaction
    }

    pZelote = zelotes[pPos] == 1;
    vZelote = zelotes[vPos] == 1;
    pIgual = polaridade[pPos] == polaridade[vPos];

    if(pZelote){
        if(!pIgual){
            zelotes[pPos] = 0;
            confianca[pPos] = 0.;
            virouN = 1;
        }
    } else {
        if(!pIgual){
            polaridade[pPos] = polaridade[vPos];
            confianca[pPos] = 0.;
            if(polaridade[pPos] == 0){
                mag--;
            } else {
                mag++;
            }
        } else {
            confianca[pPos] += deltaEta;
            if(confianca[pPos] >= 1.){
                zelotes[pPos] = 1;
                virouZ = 1;
            }
        }
    }

    // locking neighbour conviction variation
    //confianca[vPos] += deltaEta;
    if((confianca[vPos] >= 1.) && !vZelote){
        zelotes[vPos] = 1;
        if(bloqueado(vPos)){
            --nMoveis;
            update_matrix(vPos, moveis[nMoveis], quaisMoveis, moveis);
        }
        verificarLiberdadeVizinhos(vPos, 0);
    }

    if(virouZ == 1){
        if(bloqueado(pPos) && quaisMoveis[pPos] < nMoveis){
            --nMoveis;
            update_matrix(pPos, moveis[nMoveis], quaisMoveis, moveis);
        }
        verificarLiberdadeVizinhos(pPos, 0);
    }

    if(virouN == 1){
        if(quaisMoveis[pPos] >= nMoveis){
            update_matrix(pPos, moveis[nMoveis], quaisMoveis, moveis);
            ++nMoveis;
        }
        verificarLiberdadeVizinhos(pPos, 1);
    }
    return;
}

void verificarLiberdadeVizinhos(int i, int verificarLivre){
    int j, k;
    for(j = 0; j < 4; j++){
        k = vizinhos[i][j];
        if(verificarLivre == 0){
            if(bloqueado(k) && (quaisMoveis[k] < nMoveis)){
                nMoveis--;
                update_matrix(k, moveis[nMoveis], quaisMoveis, moveis);
            }
        } else {
            if(quaisMoveis[k] >= nMoveis) {
                update_matrix(k, moveis[nMoveis], quaisMoveis, moveis);
                nMoveis++;
            }
        }

    }
    return;
}



void atribuirPropriedades(void){
    int i, j = 0;
    ncl = 0;
    mag = 0;

    switch(CONFIGURACAO_I){
        case 1:break;
        case 2:
            for(i = 0; i < N; i++){
                zelotes[i] = 1;
                confianca[i] = 2.;
                if(i < N / 2){
                    polaridade[i] = 1;
                    mag++;
                } else {
                    polaridade[i] = 0;
                }
            }
            break;
        case 3:
            for (i = 0; i < N; i++){
                zelotes[i] = 1;
                confianca[i] = 2.;
                if(FRANDOM < 0.5){
                    polaridade[i] = 1;
                    mag++;
                } else {
                    polaridade[i] = 0;
                }
            }
        default: break;
    }

    for(i = 0; i < N; i++){
        if(bloqueado(i)){
            moveis[N - 1 + (j - i)] = i;//ADICIONA POLOS BLOQUEADOS NO FINAL DA LISTA DE MOVEIS
            quaisMoveis[i] = N - 1 + (j - i);
        } else {
            moveis[j] = i;
            quaisMoveis[i] = j;
            j++;
        }
    }
    nMoveis = j;
}

int bloqueado(int ii){
    int i, nPolo = polaridade[ii], nZelotes = zelotes[ii];
    for(i = 0; i < 4; i++){
        nPolo += polaridade[vizinhos[ii][i]];
        nZelotes += zelotes[vizinhos[ii][i]];
    }
    return (nPolo == 0 || nPolo == 5) && nZelotes == 5;
}


void criarVizinhos(void){
    int i,j;

    // inicializa listas
    dir = jmalloc(N * sizeof(int));
    esq = jmalloc(N * sizeof(int));
    cima = jmalloc(N * sizeof(int));
    baixo = jmalloc(N * sizeof(int));

    populacionaVizinhos();

    vizinhos = ((int **)malloc(N*sizeof(int *)));
    if(vizinhos == NULL){
        fprintf(stdout,"malloc error (neighbour)\n" );
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < N; i++){
        vizinhos[i] = jmalloc(4 * sizeof(int)); // ABRE ESPA�O PARA VIZINHOS
    }
    for(j = 0; j < N; j++){
        vizinhos[j][0] = dir[j];
        vizinhos[j][1] = esq[j];
        vizinhos[j][2] = cima[j];
        vizinhos[j][3] = baixo[j];
    }

    return;
}

void populacionaVizinhos(void){
    int i;
    for(i = 0; i < N; i++){
        if(i % L == 0) { // LADO ESQUERDO
            *(esq + i) = i + L - 1;
            *(dir + i) = i + 1;
        } else if(i % L == L - 1) { // LADO DIREITO
            *(dir + i) = i - L + 1;
            *(esq + i) = i - 1;
        } else { // MEIO
            *(dir + i) = i + 1;
            *(esq + i) = i - 1;
        }
        if(i < L) { // CIMA
            *(cima + i) = i + N - L;
            *(baixo + i ) = i + L;
        } else if(i >= N - L) { // BAIXO
            *(baixo + i ) = i - N + L;
            *(cima + i) = i - L;
        } else { // MEIO
            *(cima + i) = i - L;
            *(baixo + i ) = i + L;
        }
    }
    return;
}

void verificarCruzamentoInterface(){
    int n = 0, i, i0 = L * (L - 2) / 2, sUp, sDown, s0 = polaridade0[i0];
    if(L % 2 == 0 && FIXED_BORDERS == 1){
        for(i = 0; i < L; i++){
            sUp = polaridade0[i0 + i];
            sDown = polaridade0[i0 + i + L];
            if(sUp == s0 && sDown == s0){
                n = n + 1;
            }

        }
    }
    interfCol = n;
}


unsigned long generateFile(void){
  char name[100], context[50];
  unsigned long seed = (unsigned long) time(NULL);
  if(seed % 2 == 0) seed += 1;
  if(FORCE_SEED == 1) seed = (unsigned long) 123456789;
  switch (CONFIGURACAO_I) {
    case 1:
      sprintf(context,"PVM_L%d_IS%d_R%d",L,CONFIGURACAO_I,RAIO);
  	   break;
    case 2:
      sprintf(context,"PVM_L%d_IS%d",L,CONFIGURACAO_I);
  	  break;
    case 3:
      sprintf(context,"PVM_IS%d",CONFIGURACAO_I);
  		break;
  }
  int f = 0;
  while(f == 0){
    sprintf(name, "%s_%ld.dat", context, seed);
    fp1 = fopen(name, "r");
    if(fp1 != NULL){
      seed += 2;
    } else {
      f += 1;
    }
    fclose(fp1);
  }

  sprintf(name, "%s_%ld.dat", context, seed);
  fp1 = fopen(name, "w");

  fprintf(fp1,"# Persistent Voter Model: pvm.c\n");
  fprintf(fp1,"# Initial state: %d\n",CONFIGURACAO_I);
  fprintf(fp1,"# Linear size: %d",L);
  if (CONFIGURACAO_I==1) fprintf(fp1,"  Radius: %d\n",RAIO);
                   else fprintf(fp1,"\n");
  fprintf(fp1,"# Maximum simulation time: %d\n",(int)TEMPO_MAX);
  fflush(fp1);

  return seed;
}

void save_configuration(int qual, char* name){
int i;
char output_file_eps[100];

sprintf(output_file_eps,"PVM_%s_%06d.eps",name, qual);
lat2eps_init(L, L);
//lat2eps_set_color(0,0x000000); // black, for s=0, normal sites
//lat2eps_set_color(1,0x31a2f2); // blue, for s=0, zealot sites
//lat2eps_set_color(2,0xFFFFFF); // white, for s=1, normal sites
//lat2eps_set_color(3,0xFF0000); // red, for s=1, zealot sites
lat2eps_set_color(0,0x0000FF); // dark blue, for s=0, normal
lat2eps_set_color(1,0x81C2EF); // light blue, for s=0, zealot
lat2eps_set_color(2,0xFF0000); // dark red, for s=1, normal
lat2eps_set_color(3,0xFA8282); // light red, for s=1, zealot
for (i=0; i<N; ++i)
    lat2eps_set_site(i%L,i/L,2*polaridade0[i]+zelotes0[i]);
lat2eps_gen_eps(output_file_eps,0,0,L,L,1,3);
lat2eps_release();

return;
}

void unique_cluster(){
  int *label, *z0c, *z1c, *arc, *normal_label;
  int sV1, sV2, zV1, zV2, s0, z0;
  label = jmalloc(N * sizeof(int));
  normal_label = jmalloc(N * sizeof(int));
  z0c = jmalloc(N * sizeof(int));
  z1c = jmalloc(N * sizeof(int));
  arc = jmalloc(N * sizeof(int));
  for (int i = 0; i < N; ++i){
      label[i] = i;           //label are from cluster with same zealotery and spins
      normal_label[i] = i;    //normal_label are from clusters with same zealotery
  }


  //diferencia os clusters
  for (int i = 0; i < N; ++i) {
    s0 = polaridade0[i]; z0 = zelotes0[i];
    sV1 = polaridade0[dir[i]]; zV1 = zelotes0[dir[i]];
    sV2 = polaridade0[baixo[i]]; zV2 = zelotes0[baixo[i]];
    if(z0 == zV1){
      unionfind(i, dir[i], normal_label);
      if(s0 == sV1){
        unionfind(i, dir[i], label);
      }
    }
    if(z0 == zV2){
      unionfind(i, baixo[i], normal_label);
      if(s0 == sV2){
        unionfind(i, baixo[i], label);
      }
    }
  }
  if(L == 24){
    for(int i = 0; i < L; i++){
      for(int j = 0; j < L; j++){
        int index = (i * L) + j;
        printf("%d   ", label[index]);
      }
      printf("\n");
    }
  }

  int z0_lab = label[0], z1_lab = label[N - 1], s = 0;
  for(int i = 0; i < N; i++){
    if(z0_lab == label[i]){
      z0c[i] = 1;
    }
    if(z1_lab == label[i]){
      z1c[i] = 1;
      if(L == 16){
        printf("%d\n", i);
        s += 1;
      }
    }
  }

  getBiggestCluster(arc, normal_label, 2);

  l0_int = 0;
  l1_int = 0;
  l0_int = getInterfacialSizeOfCluster(z0c, arc);
  l1_int = getInterfacialSizeOfCluster(z1c, arc);
  free(label);free(z0c); free(z1c);
}

int getInterfacialSizeOfCluster(int *cluster1, int *cluster2){
  int l = 0, si, zi, viz, lab;
  for(int i = 0; i < N; i++){
    lab = cluster1[i];
    if(lab != 1) continue;
    si = polaridade0[i]; zi = zelotes0[i];

    for(int j = 0; j < 4; j++){
      viz = vizinhos[i][j];

      if(polaridade0[viz] != si || zelotes0[viz] != zi){
        l += 1;
      }


    }
  }
  return l;

}

void getBiggestCluster(int *biggestCluster, int *label, int config){
  int *clusterSize_arr = jmalloc(N * sizeof(int));
  int j;
  for(int i = 0; i < N; i++) clusterSize_arr[i] = 0;
  for(int i = 0; i < N; i++){
    j = label[i];
    if(config > 1){
      if(zelotes0[i] == 0){
        clusterSize_arr[j] += 1;
      }
    } else {
      if(polaridade0[i] == config){
        clusterSize_arr[j] += 1;
      }
    }
  }
  int clusterSize = 0, cluster_lab = 0;
  for(int i = 0; i < N; i++){
    if(clusterSize <= clusterSize_arr[i]){
      clusterSize = clusterSize_arr[i];
      cluster_lab = i;
    }
  }
  for(int i = 0; i < N; i++){
    if(label[i] == cluster_lab){
      biggestCluster[i] = 1;
    } else {
      biggestCluster[i] = 0;
    }
  }
  free(clusterSize_arr);

}


void hoshen_kopelman(int *his0, int *numc, int *per){
    int i,j,*label,*siz;
    int temp1,temp2,wtemp1,wtemp2,probpercv0,probpercv1;

    //criar uma lista label de 0 a N identificando clusters
    //a lista siz guarda o tamanho da cluster i-esima
    label = jmalloc(N * sizeof(int));
    siz = jmalloc(N * sizeof(int));
    for (i = 0; i < N; ++i){
        label[i] = i;
        siz[i] = 0;
    }

    probpercv0=0;
    probpercv1=0;

    //diferencia os clusters
    for (i = 0; i < N; ++i) {
        if ( (2 * polaridade0[i] + zelotes[i]) == (2 * polaridade0[dir[i]] + zelotes0[dir[i]])) unionfind(i, dir[i], label);
        if ( (2 * polaridade0[i] + zelotes[i]) == (2 * polaridade0[baixo[i]] + zelotes0[baixo[i]])) unionfind(i, baixo[i], label);
    }

    for (i = 0; i < N; ++i){                     /* Measure the cluster sizes */
         j = i;
         while (label[j] != j)                  /* it doesn't point to itself */
               j = label[j];                    /* follow the cluster */
         ++siz[label[j]];                       /* update the cluster size */
         label[i] = label[j];                   /* update the cluster label */
    }

    temp1 = -1;                                /* Size of the biggest cluster */
    temp2 = -1;                                /* Size of the second biggest cluster */
    wtemp1 = 0;
    wtemp2 = 0;
    *numc = 0;                                  /* number of clusters */

    for (i = 0; i < N; ++i) {
        if (siz[i] > 0){
            if (siz[i] >= L){
                if (polaridade0[i] == 0) probpercv0 += percolates2d(i, label);
                if (polaridade0[i] == 1) probpercv1 += percolates2d(i, label);
            }
            if (siz[i] >= temp1) {
                temp2 = temp1;
                temp1 = siz[i];
                wtemp2 = wtemp1;
                wtemp1 = i;
            } else if (siz[i] > temp2) {
                temp2 = siz[i];
                wtemp2 = i;
            }
            ++(*numc); /* Count total number of clusters */
          }
    }

    *per = probpercv0 + probpercv1;
    zar = 0;
    for (i = 0; i < N; ++i)
        if ((i != wtemp1) && (i != wtemp2) && (zelotes[i]==1))
           {
            zar += siz[i];

           }
    //zar MEDE QUANTOS ZELOTES EXISTEM NA AR, comparar agora com o tamanho m�dio dos clusters na AR (sem separar os zelotes)... a AR deve parar de crescer quando o raio m�dio desses clusters fica maior do que uma parede viaja em 1/deta passos

    free(siz); /* from here on, we use the smaller newsiz */
    free(label);
    return;
}



int percolates2d(int qual,int *lab) {
    int ok1=0,ok2=0;
    int i,j;
    for (i = 0; i < N; i += L){
        ok2 = 0;
        for (j = i; j < i + L; ++j){
            if (lab[j] == qual) {
                ok2=1;
                break;
            }
        }
        if (ok2 == 0) break;
    }
    return ok1 + ok2;
}
