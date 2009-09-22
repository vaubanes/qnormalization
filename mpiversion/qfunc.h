/**
* qfunc.h
*
*  Author : José Manuel Mateos
*
**/

#ifndef QFUNC_H_
#define QFUNC_H_

#define BLOCK_SIZE 		 		2
#define MAX_SIZE_LINE   		500
#define DEFAULT_NUM_GENES       15
#define DEFAULT_NUM_EXPERIMENTS 2

#ifdef DEBUG
  #define info(s) fprintf(stderr,"INFO: %s\n",s);
#else
  #define info(s) 
#endif

typedef struct { // List of files
  char *FileName;
  int NumGenes;
  char FileType;
} InfoFile;

typedef struct { // Contains Average of array
  double Value;
  int Elements;
} Average;

typedef struct {
	int  BlockSize;                  	  		// Size of block
	char FileListExperiments[MAX_SIZE_LINE];    // File with a list of experiments files
	char FileOut[MAX_SIZE_LINE];    			// Output file name
	int  Traspose;            		     		// Traspose file to file final results (0:NOT 1:Yes)
	int  NumGenes;                       		// Number of Genes (rows)
	int  NumExperiments;                  		// Number of Experiments or samples(cols)
} Params;


// Function prototypes------------------------------------------

// general functions
Params *commandline(int argc, char *argv[]);
InfoFile* load_input_files(Params *);
void load_parcial_result(InfoFile*flist, int col, double *dataIn, int num_genes);
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

