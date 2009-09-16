/**
* qfunc.h
*
*  Author : José Manuel Mateos
*
**/

#ifndef QFUNC_H_
#define QFUNC_H_

#define BLOCK_SIZE 2
#define MAXLIN   500
#define NGEN     15
#define NEXP     2

#ifdef DEBUG
  #define info(s) fprintf(stderr,"INFO: %s\n",s);
#else
  #define info(s) 
#endif

typedef struct  { // LIst of files
  char *fname;
  int num_genes;
  char fType;
}Files;

typedef struct { // Average array
  double av;
  int num;
} Average;

typedef struct {
	int  block_size;                  // Size of block
	char flist_experiments[MAXLIN];   // file with a list of files
	char fout[MAXLIN];    			  // Output file name
	int  traspose;            		  // Traspose file to file final results (0:NOT 1:Yes)
	int  n_genes;                     // Number of Genes (rows)
	int  n_experiments;               // Number of Experiments or samples(cols)
} Params;


// Function protorypes------------------------------------------

// general functions
Params *commandline(int argc, char *argv[]);
Files* load_input_files(Params *);
void load_parcial_result(Files*flist, int col, double *dataIn, int num_genes);
void terror(char *);


void debug_print(char *, double*, int);
int  transpose_matrix(Params*);
void quick_sort(double *array,int l,int r,int *index);
int partition( double* a, int l, int r, int *indexes);


// related to Qnorm

int accumulate_row(Average *average, double *input , int num_genes);
int qnorm(double *input, int *dindex, int num_genes);
int calculate_blocks(int num_experiments, int block_size);
int calculate_initial_blocks(int numBlocks, int nProcesors);
int calculate_index_blocks(int *index,int num_experiments,int block_size);

#endif /* QFUNC_H_ */

// ===============================================================================

