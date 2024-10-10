/********************************************************************
***                       Monte Carlo Routines                    ***
***                         V0.01 17/01/2024                      ***
***                                                               ***
***                                                               ***
*** List of functions (auxiliary routines, transparent to the     ***
*** user, are not listed:                                         ***
***                                                               ***
*** FRANDOM: return a pseudorandom number in [0,1)                ***
*** start_randomic: returns a seed for the random generator,      ***
***                 either a time or a fixed value (for           ***
***                 debugging purposes)                           ***
*** ngaussian: gaussian random values (CHECK)                     ***
*** jmalloc:                                                      ***
*** neighbours_int:                                               ***
*** update_matrix:                                                ***
*** unionfind:                                                    ***
*** create_time_table:                                            ***
*** initial_configuration_plusminus: (void)                       ***
*** sum_1d_array: sum all elements of a 1d array (int)            ***
*** hamming_distance: (int)                                       ***
*** ***
*** ***
*** ***
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> 
#include <string.h>
#include <time.h>  

#define MC_VERSION  "0.01"
#define FNORM   (2.3283064365e-10)
#define RANDOM  ( (mc_ira[ip++]=mc_ira[ip1++]+mc_ira[ip2++]) ^mc_ira[ip3++] )
#define FRANDOM (FNORM * RANDOM)

/********************************************************************
***                      Variable Declarations                    ***
********************************************************************/

unsigned mc_zseed, mc_ira[256];
unsigned char ip,ip1,ip2,ip3;

/********************************************************************
*            Random Number Generator by Parisi & Rapuano            *
*                  Last Modified: 18/05/2000                        *
*                                                                   *
*  First, the function start_randomic() should be called to create  *
*  the seed (odd number):                                           *
*                 seed = start_randomic();                          *
*                                                                   *
*  Return: start_randomic -> unsigned long int                      *
*          FRANDOM -> double in [0,1)                               *
********************************************************************/
unsigned rand4init(void)
{
     unsigned long long y;

     y = (mc_zseed*16807LL);
     mc_zseed = (y&0x7fffffff) + (y>>31);
     if (mc_zseed&0x80000000)
         mc_zseed = (mc_zseed&0x7fffffff) + 1;
     return mc_zseed;
}

void Init_Random(void)
{
     int i;

     ip=128;
     ip1=ip-24;
     ip2=ip-55;
     ip3=ip-61;

     for (i=ip3; i<ip; i++)
         mc_ira[i] = rand4init();
}

/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 13/08/2018                    *
*                                                               *
* When debugging, it uses always the same seed, otherwise,      *
* it is the number of seconds ellapsed since 00:00:00 UTC       *
* 01/01/1970.                                                   *
****************************************************************/
void start_randomic(long semente)

{
if (semente==0) semente = (long) time(NULL);      /* random seed */
if (semente%2==0) ++semente; /* odd number */
#ifdef DEBUG
 semente = 123456789;
#endif
mc_zseed = semente;
Init_Random();
return;
}

/********************************************************************
***                   Gaussian Random Number Generator            ***
***                     Last Modified: 18/05/2000                 ***
********************************************************************/ 
float ngaussian(void)

{
static int iset=0;
static float gset;
float fac,r,v1,v2;
   
if (iset==0) {
	      do {
	  	  v1=2.0*FRANDOM-1.0;
		  v2=2.0*FRANDOM-1.0;
		  r=v1*v1+v2*v2;
	         } 
	      while (r>=1 || r==0.0);
	      fac=sqrt(-2.0*log(r)/r);
	      gset=v1*fac;
	      iset=1;
	      return v2*fac;
             }
        else {
	      iset=0;
	      return gset;
             }
}

/**********************************************************************************
 *                                   malloc                                       *
 *                        Last modified: 19/07/2006                               *
 *********************************************************************************/
void *jmalloc (unsigned int nbytes)
  
{
   
void *ptr;
ptr = malloc (nbytes);
if (ptr == NULL)
   {
    printf("malloc returned NULL!\n");
    exit(EXIT_FAILURE);
   }
return ptr;
}

/*****************************************************************
***                          3D Neighbours                     ***
***                    Last Modified: 10/01/2024               ***
***                                                            ***
***  In a 3D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
***  Input: right,left,up,down,front, back -> neighb. sites    ***
***         lsize -> linear dimension (lsize^3)                ***
***         dimension -> 1, 2 or 3                             ***
***                                                            ***
*** For d<3, some of the above matrices may be NULL.           ***
***                                                            ***
*****************************************************************/
void neighbours_int(int *right,int *left,int *up, int *down,int *front, 
                    int *back,int lsize, int dimension)

{
int i,l2,l3=lsize;

l2 = lsize*lsize;
if (dimension>1) l3 *= lsize;
if (dimension==3) l3 *= lsize; 

for (i=0; i<l3; ++i) 
    {
     if (i % lsize==lsize-1) *(right+i) = i - lsize + 1; /* righmost plane */
                        else *(right+i) = i+1;
     if (i % lsize==0) *(left+i) = i + lsize - 1;        /* leftmost plane */
                  else *(left+i) = i-1;
     if (dimension>1)
        {
         if (i % l2 < lsize) *(up+i) = l2 - lsize + i;       /* top plane */
                        else *(up+i) = i - lsize;
	 if (i % l2 >= (l2-lsize)) *(down+i) = i + lsize - l2;  /* bottom plane */
	                      else *(down+i) = i + lsize;
	}
     if (dimension==3)
        {
         if (i < l2) *(front+i) = l3 - l2 + i;
                else *(front+i) = i - l2; /* frontal plane */
         if (i >= (l3-l2)) *(back+i) = i % l2;
                       else *(back+i) = i + l2; /* back plane */
	}
    }
return;
}         


/*********************************************************************
***                     Update Localization Matrixes               ***
*********************************************************************/
void update_matrix(int s1, int s2, int *which_emp, int *emp)

{
int sitetmp;

emp[which_emp[s1]] = s2;
emp[which_emp[s2]] = s1;
sitetmp = which_emp[s1];
which_emp[s1] = which_emp[s2];
which_emp[s2] = sitetmp;
return;
}

/*******************************************************************************
*                   Instantaneous Percolation Measurements                     *
*                     Last modified:  17/01/2024                               *
*                                                                              *
* Remember that a cluster is identified by the last site in the chain of       *
* pointers.                                                                    *
*******************************************************************************/
void unionfind(int i,int j, int *lab)
          
{
int i1,j1;

i1 = lab[i];                            /* check where i points to          */
while (lab[i1] != i1)                   /* while it doesn't point to itself */
      i1 = lab[i1];

j1 = lab[j];                            /* check where j points to          */
while (lab[j1] != j1)                   /* while it doesn't point to itself */
      j1 = lab[j1];

if (lab[i1] > lab[j1]) lab[i1] = lab[j1];
                  else lab[j1] = lab[i1];
return;
}

/**********************************************************************
***                            Time Table                           ***
***                    Last Modified: 16/01/2024                    ***
***  The total number of time steps and the number of measures are  ***
***  specified. The last argument tells whether the scale is linear ***
***  or logarithmic.                                                ***
**********************************************************************/
void create_time_table(int *t1,int total_time, int measures, int linear_or_log)

{
int i,temp0;
double temp,temp1;

switch (linear_or_log)
       {
        case 0: // linear intervals of time
                t1[0] = 0; /* thermalization time or t=0*/
                temp0 = total_time/measures;
                for (i=1; i<=measures; ++i)
                    t1[i] = t1[i-1] + temp0;
                break;
        case 1: // logarithmic intervals         
                t1[0] = 0;
                t1[1] = 1;
                temp = (double) t1[1];  
                temp1 = pow(((double)total_time)/t1[1],1./(measures-1));                
                for (i=2; i<=measures; ++i)
                    {
                     temp *= temp1;
                     t1[i] = (int) temp;
                    }
                break;
        default: printf("The last argument of create_time_table may be either 0 or 1!\n");
                 exit(1);
       }
return;
}

/************************************************************************
***                     Create Initial Configuration                  ***
***                      Last modified: 16/01/2024                    ***
***                                                                   ***
*** Create a random initial configuration with a given magnetization. ***
*************************************************************************/
void initial_configuration_plusminus(float mag, int ll, int *s)

{
int i;
float rm;
   
rm = 0.5*(1.0 - mag);  /* fraction of -1 sites */
for (i=0; i<ll; ++i)
    {
     if (FRANDOM<rm) s[i] = -1;
                else s[i] =  1;
    }
return;
}

/*********************************************************************
***                          Sum 1d array                          ***
***                   Last Modified: 17/01/2024                    ***
***                                                                ***
*** Return: LSIZE*magnetization                                    ***
*********************************************************************/
int sum_1d_array(int *s, int lsize)

{
unsigned long i,temp=0;
   
for (i=0; i<lsize; ++i)
    temp += s[i];
return temp;
}

/**********************************************************************
***                         Hamming Distance                        ***
***                   Last modified: 24/01/1999                     ***
***    Returns the number of different sites between two vectors.   ***
**********************************************************************/ 
int hamming_distance(int *s1, int *s2, int lsize)

{
unsigned long i,temp=0;
  
for (i=0; i<lsize; ++i)
     if (s1[i] != s2[i]) ++temp;
return temp;
}
