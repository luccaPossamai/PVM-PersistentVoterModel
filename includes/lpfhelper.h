#include <stdio.h>
#include <stdlib.h>

/****************LPFileHelper*****************
    Make writable file more easy to manage
*********************************************/


FILE* safeOpen(char*, char*);
FILE* safeSeedOpen(char*, char*, unsigned int*, int);
int existFile(char*);

/*****************safeOpen()******************
    Open a file with 'name'_0.'extension'
    Saves with an int count at the end. 
*********************************************/
FILE* safeOpen(char *name, char *extension){
	FILE *f;
	char name0[50];
	int i = 0;
	do {
		sprintf(name0, "%s_%d%s",name, i, extension);
		i += 1;
	} while (existFile(name0) == 1);
	f = fopen(name0, "w");
	return f;
}

/****************existFile()*****************
    Verify if the file with name exists
    Require extension
*********************************************/
int existFile(char *name){
	FILE *f;
	int i = 0;
	f = fopen(name, "r");
	if(f != NULL) {
		i = 1;
		fclose(f);
	}
	return i;
}

/***************safeSeedOpen()****************
    Returns a file with an odd seed;
    If the file exists return the next free 
    odd seed file;
*********************************************/
FILE* safeSeedOpen(char *name, char *extension, unsigned int *seed, int forceSeed){
	FILE *f;
	char name0[150];
	if(*seed % 2 == 0) *seed += 1;
	do {
		sprintf(name0, "%s_S%d%s",name, *seed, extension);
		*seed += 2;
		
	} while (existFile(name0) == 1 && forceSeed == 0);
	f = fopen(name0, "w");
	return f;
}


