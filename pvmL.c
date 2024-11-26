#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mc2023.h>
//#include <lat2eps.h>
#include <monte_carlo.h>
#include <fhelper.h>
#include <stdbool.h>


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
unsigned int generateFile(void);
void hoshen_kopelman(int*, int*, int*);
void verificarLiberdadeVizinhos(int, int);
bool percolates2d(int, int*);
void createPerfectLogTable(int*);
void comecar(int);
void iniciar_coleta(int);
void verificarCruzamentoInterface(void);
void save_configuration(int, char*);
void time_evolution(void);
void write(double, int, int, int);
void updateOldValues(void);
void interface(void);
int getInterfacialSizeOfCluster(int*);
void getBiggestCluster(int*, int*, int);
void unique_cluster(void);
void clusterNumber(int*, int*);


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
#define MEDIDAS                 20

/********** Par�metros para EVOLUCAO = 1(consensus) ****************/
#define DELTAETA_FINAL          1.0
#define DELTAETA_MEDIDAS        41

/******* Par�metros para CONFIGURACAO_I = CIRCULAR *****************/
#define RAIO                    10          //NAO IMPLEMENTADO AINDA

#define FORCE_SEED              1
/******* Par�metros para CONDICAO DE CONTORNO LINEAR FIXA **********/
#define FIXED_BORDERS			1			//non-periodic borders don't interact
#define B						1			//"b" = "boundry condition"  
											//0 -> periodic
											//1 -> dobrushin vertical(periodic horizontally)
											//2 -> dobrushin horizontal(periodic vertically)
											//3 -> square		   
#define L_INICIAL               16
#define L_FINAL                 256
/*******************************************************************/
#define SAVE_CONFIG             1
#define CONFIG_T0               0	
#define CONFIG_DELTA_T          100000
#define CONFIG_TF               111000//CONFIG_T0 + 10 * CONFIG_DELTA_T



int *cima, *baixo, *esq, *dir, **vizinhos, *zelotes, *polaridade/*0, 1 */,
    *moveis, *quaisMoveis, *tempos, mag, ncl, zar, interfCol;
int nMoveis0, mag0, *polaridade0, *zelotes0,  *moveis0, *quaisMoveis0,
    l0_int = 0, l1_int = 0, L = LSIZE, N = NSIZE;
float *confianca0;
int nMoveis;
double tempo,deltaEta;
float *confianca, *o_confianca;
unsigned int seed;
FILE *fp1;

int main(void){
    int lps = 1;//getIntInputFromClient("loops = ", 1, 1e4);
    if(FORCE_SEED == 0){
        lps = getIntFromUser("Measures", 1, 1e4);
    }
    iniciar_coleta(lps);

}

void iniciar_coleta(int i){
    int j, t = 0;
    seed = 0;
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
	seed = generateFile();
	setupRandom(seed);
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
    fclose(fp1);
}

void teste(void){

}



void post(void){
    free(dir);
    free(esq);
    free(cima);
    free(baixo);
    free(vizinhos);
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

}



void pre(){
	criarPropriedades();
	createNeighboursMatrix(vizinhos, N);

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
    vizinhos = smalloc(N * sizeof(int*));

}

void updateOldValues(void){
    nMoveis0 = nMoveis;
    mag0 = mag;
    for(int i = 0; i < N; i++){
    	zelotes0[i] = zelotes[i];
		confianca0[i] = confianca[i];
		polaridade0[i] = polaridade[i];
		moveis0[i] = moveis[i];
		quaisMoveis0[i] = quaisMoveis[i];
    }
    
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
	float *l_arr = jmalloc(MEDIDAS * sizeof(float));
	geomProgression(l_arr, L_INICIAL, L_FINAL, MEDIDAS);

	char name[50];
	
    int threshold = 1e3, l0, l1, lar;
    deltaEta = DELTAETA_INICIAL;
    threshold = 1e3;
    fprintf(fp1,"#  deltaEta: %.5f\n", deltaEta);
    fprintf(fp1,"#  threshold: %d\n", threshold);
    fprintf(fp1,"#  L    <L0>   <L1>   <L_Ar>\n");
    for(int i = 0; i < MEDIDAS; i++){
      L = (int) l_arr[i];
      N = L * L;
        sprintf(name, "%d", L);
        post();
        pre();
        atribuirPropriedades();
        updateOldValues();
        tempo = 0.0;
        l0 = 0; l1 = 0; lar = 0;
        while(tempo < TEMPO_MAX){
          tempo += 1./nMoveis;
            interacoes();
            if(tempo >= threshold){
              unique_cluster();
              l0 = l0_int;
              l1 = l1_int;
              fprintf(fp1,"%d %d %d %d\n", L, l0 - L, l1 - L, lar);
              fflush(fp1);

              if(SAVE_CONFIG) save_configuration(tempo, name);
              break;
            }
            if(tempo + 1./nMoveis >= threshold){
            	updateOldValues();
            }
            if(tempo > TEMPO_MAX){
                break;
            }
        }
    }
}

void write(double proxTempo, int z, int z1, int per){
    fprintf(fp1,"%.5f  %d %d %d %d %d %d %d %d\n", proxTempo, mag0, z1, z, ncl, zar, per, nMoveis0, interfCol);
    if(SAVE_CONFIG == 1){
        if (proxTempo >= CONFIG_T0 && proxTempo <= CONFIG_TF) {
            //save_configuration(proxTempo, "evolution");
        }
 
    }
    fflush(fp1);
}

void time_evolution(void){
    double proxTempo = 0;
    float *time_arr;
    int z = 0, z1 = 1, per;
    time_arr = smalloc(MEDIDAS * sizeof(float));
    if(ESCALA == 0){
    	linearProgression(time_arr, 0, TEMPO_MAX, MEDIDAS);
    } else {
	    geomProgression(time_arr, 0, TEMPO_MAX, MEDIDAS);    
    }

    deltaEta = DELTAETA_INICIAL;
    fprintf(fp1,"# Deta = %.5f\n",deltaEta);
    fprintf(fp1,"#  t  m  z1  z  ncl z_ar per mob int_cross\n");
    atribuirPropriedades();
    updateOldValues();
    tempo = 0.0;
	for(int i = 0; i < MEDIDAS; i++){
		proxTempo = time_arr[i];
		while(nMoveis > 0 && tempo < proxTempo){
			tempo += 1./nMoveis;
			interacoes();
			if(tempo + 1./nMoveis >= proxTempo){
				updateOldValues();
			}
		}
		z = 0;
        z1 = 0;
        for(int i = 0; i < N; i++){
        	if(polaridade0[i] == 1){
        		z1 += zelotes0[i];
            }
            z += zelotes0[i];
        }
        clusterNumber(&ncl,&per);
        verificarCruzamentoInterface();
		write(proxTempo, z, z1, per);
	}
    return;
}

void interacoes(void){
    int pPos,pIndice, vIndice, vPos, virouZ = 0, virouN = 0, pZelote, vZelote, pIgual;

    pIndice = randomInt(nMoveis);
    pPos = moveis[pIndice];
    vIndice = randomInt(4);
    vPos = vizinhos[pPos][vIndice];

	if(!isValidInteraction(pPos, vPos, N, B)) return;
	if(FIXED_BORDERS == 1 && isOnNonPeriodicBorder(pPos, N, B)) return;
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

void verificarCruzamentoInterface(){
    int n = 0, i, i0 = L * (L - 2) / 2, sUp, sDown, s0 = polaridade0[i0];
    if(L % 2 == 0){
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


unsigned int generateFile(void){
  char name[200], context[100];
  unsigned int id = setupRandom(0);
  if(FORCE_SEED == 1) id = (unsigned int) 123456789;
  sprintf(context, "data_pvm_C%d_BC_%d", EVOLUCAO, B);
  switch(EVOLUCAO){
  	case 0:
  		sprintf(name, "%s_L%d_I%d", context, L, CONFIGURACAO_I);
  	case 1:
  		sprintf(name, "%s_L%d", context, L);
  	case 2:
  		sprintf(name, "%s_dE%.4f", context, deltaEta);
  }
  fp1 = safeSeedOpen(context, ".dat", &id, FORCE_SEED);
  fprintf(fp1,"# Persistent Voter Model: pvm.c\n");
  //fprintf(fp1,"# Initial state: %d\n",CONFIGURACAO_I);
  //fprintf(fp1,"# Linear size: %d",L);
  //if (CONFIGURACAO_I==1) fprintf(fp1,"  Radius: %d\n",RAIO);
    //               else fprintf(fp1,"\n");
  //fprintf(fp1,"# Maximum simulation time: %d\n",(int)TEMPO_MAX);
  fflush(fp1);

  return id;
}

void save_configuration(int qual, char* name){
/*
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
*/
}

void unique_cluster(){
  int *label, *lat, *z0c, *z1c;
	  
	label = jmalloc(N * sizeof(int));
	lat = smalloc(N * sizeof(int));
	z0c = smalloc(N * sizeof(int));
	z1c = smalloc(N * sizeof(int));
	for(int i = 0; i < N; i++){
		lat[i] = polaridade0[i] + 2 * zelotes0[i];
		z0c[i] = -1;
		z1c[i] = -1;
	}
	  
  	labelCluster(label, lat, vizinhos, N, 1);
	int z0_lab = label[0], z1_lab = label[N - 1];
  	for(int i = 0; i < N; i++){
		if(z0_lab == label[i]){
		  z0c[i] = 1;
		}
		if(z1_lab == label[i]){
		  z1c[i] = 1;
		}
  	}
  	l0_int = 0;
  	l1_int = 0;
  	l0_int = getInterfacialSizeOfCluster(z0c);
  	l1_int = getInterfacialSizeOfCluster(z1c);
  	if(L == 16){
  		for(int i = 0; i < L; i++){
  			for(int j = 0; j < L; j++){
  				printf("%d ", lat[i * L + j]);
  			}
  			printf("\n");
  		}
  	}
  	free(label);free(z0c); free(z1c);
}

int getInterfacialSizeOfCluster(int *cluster1){
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

void clusterNumber(int *nC, int *nPer){
	int *label, *size, *lat;
	lat = smalloc(N * sizeof(int));
	label = smalloc(N * sizeof(int));
	size = smalloc(N * sizeof(int));
	for(int i = 0; i < N; i++){
		lat[i] = polaridade0[i] + 2 * zelotes0[i];
		size[i] = 0;
	}
	labelCluster(label, polaridade0, vizinhos, N, B);
	for(int i = 0; i < N; i++){
		size[label[i]] += 1;
	}
	*nC = 0;
	*nPer = 0;
	for(int i = 0; i < N; i++){
		if(size[i] > 0){
			*nC += 1;
			if(percolates2d(label[i], label)) *nPer += 1;
		}
	}
	free(label);
	free(lat);
	free(size);	
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
    /*
    for (i = 0; i < N; ++i) {
        if ( (2 * polaridade0[i] + zelotes[i]) == (2 * polaridade0[dir[i]] + zelotes0[dir[i]])) unionfind(i, dir[i], label);
        if ( (2 * polaridade0[i] + zelotes[i]) == (2 * polaridade0[baixo[i]] + zelotes0[baixo[i]])) unionfind(i, baixo[i], label);
    }

*/
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



bool percolates2d(int qual,int *lab) {
    bool ok1 = false,ok2 = false;
    for(int i = L; i < N; i += L){
    	if(lab[i - 1] == qual) ok1 = true;
    	if(lab[i] == qual) ok2 = true;
    	if(ok1 && ok2) break;
    }
    return ok1 && ok2;
}
