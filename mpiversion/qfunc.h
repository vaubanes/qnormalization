/******************************************************************************
* FILE: Qfunc.h
* ots@ac.uma.es
* LAST REVISED: 23/02/09
******************************************************************************/

#pragma once

#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <math.h>

#define BLOCK_SIZE 2
#define MAXLIN   500
#define MAXnP    4
#define NGEN     15
#define NEXP     2
#define max(a,b)    (((a)>(b)) ? (a):(b))
#define min(x,y)    (((x) < (y)) ? (x) : (y))


struct Files { // LIst of files
  char *fname;
  int nG;
  char fType;
};

struct Average { // Average array
  double Av;
  int num;
};


struct params {// Parameters struct-----------------
  int  nP;                  // number of nodes
  char fListName[MAXLIN];   // file with a list of files
  char fOutName[MAXLIN];    // Output file name
  int  Traspose;            // Traspose file to file final results (0:NOT 1:Yes)
  int  MemIndex;            // store Index in (1) memory or (0) in disk
  int  nG;                  // Number of Genes (rows)
  int  nE;                  // Number of Experiments or samples(cols)
};

// Function protorypes------------------------------------------

// general functions
struct params *CommandLine(int argc, char *argv[]);
struct Files* LoadListOfFiles(struct params *);
void LoadFile(struct Files*, int, double *,int nG);
void terror(char *);
void Alerta(char *,char *);

void DebugPrint(char *, double*, int);
int  TransposeBin2Txt(struct params*);
void QsortC(double *array,int l,int r,int *index);
int partition( double* a, int l, int r, int *indexes);

// related to Qnorm
void QNormMain(struct params*, struct Files*);
int AccumulateRow(struct Average *, double *, int);
int Qnorm1(double *, int *, int);
int calculateBlocks(int nE);
int calculateInitialBlocks(int numBlocks, int nProcesors);
int calculateIndexBlocks(int *index,int nE);



// ===============================================================================

