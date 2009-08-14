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


struct Files { // LIst of files
  char *fname;
  int num_genes;
  char fType;
};

struct Average { // Average array
  double av;
  int num;
};


struct Params {// Parameters struct-----------------
  int  block_size;                  // Size of block
  char flist_experiments[MAXLIN];   // file with a list of files
  char fout[MAXLIN];    			// Output file name
  int  traspose;            		// Traspose file to file final results (0:NOT 1:Yes)
  int  n_genes;                  	// Number of Genes (rows)
  int  n_experiments;               // Number of Experiments or samples(cols)
};

// Function protorypes------------------------------------------

// general functions
struct Params *commandline(int argc, char *argv[]);
struct Files* load_input_files(struct Params *);
void load_parcial_result(struct Files*flist, int col, double *dataIn, int num_genes);
void terror(char *);


void debug_print(char *, double*, int);
int  transpose_matrix(struct Params*);
void quick_sort(double *array,int l,int r,int *index);
int partition( double* a, int l, int r, int *indexes);


// related to Qnorm

int accumulate_row(struct Average *average, double *input , int num_genes);
int qnorm(double *input, int *dindex, int num_genes);
int calculate_blocks(int num_experiments, int block_size);
int calculate_initial_blocks(int numBlocks, int nProcesors);
int calculate_index_blocks(int *index,int num_experiments,int block_size);

#endif /* QFUNC_H_ */

// ===============================================================================

