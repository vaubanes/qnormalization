/**
* qfunc.h
*
*  Author : José Manuel Mateos
*
**/

#ifndef QFUNC_H_
#define QFUNC_H_



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
  int num_genes;
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
struct params *commandline(int argc, char *argv[]);
struct Files* load_input_files(struct params *);
void load_parcial_result(struct Files*flist, int col, double *dataIn, int num_genes);
void terror(char *);


void debug_print(char *, double*, int);
int  transpose_matrix(struct params*);
void sort(double *array,int l,int r,int *index);
int partition( double* a, int l, int r, int *indexes);

// related to Qnorm

int accumulate_row(struct Average *average, double *input , int num_genes);
int qnorm(double *input, int *dindex, int num_genes);
int calculate_blocks(int nE);
int calculate_initial_blocks(int numBlocks, int nProcesors);
int calculate_index_blocks(int *index,int num_experiments);

#endif /* QFUNC_H_ */

// ===============================================================================

