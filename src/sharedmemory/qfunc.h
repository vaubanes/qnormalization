/******************************************************************************
* FILE: Qnorm.h
* ots@ac.uma.es
* LAST REVISED: 23/02/09
******************************************************************************/


#ifndef QFUNC_H_
#define QFUNC_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_SIZE_LINE   1000
#define DEFAULT_PROCESSOR_NUMBER    4
#define DEFAULT_NUM_GENES     15
#define DEFAULT_NUM_EXPERIMENTS     2
#define max(a,b)    (((a)>(b)) ? (a):(b))
#define min(x,y)    (((x) < (y)) ? (x) : (y))


struct files { // LIst of files
  char *fname;
  int num_genes;
  char fType;
  int pos;
};

struct average { // Average array
  float average;
  int num;
};


struct params {// Parameters struct-----------------
  int  num_processors;      // number of nodes
  char file_list[MAX_SIZE_LINE];   // file with a list of files
  char file_out[MAX_SIZE_LINE];    // Output file name
  int  transpose;            // Traspose file to file final results (0:NOT 1:Yes)
  int  mem_index;           // store Index in (1) memory or (0) in disk
  int  num_genes;           // Number of Genes (rows)
  int  num_experiments;     // Number of Experiments or samples(cols)
  int  verbose;             // Not(0) (default) / yes (1)
};

// Function protorypes------------------------------------------

// general functions
struct params *command_line(int, char **);
struct files* load_list_files(struct params *);
void load_file(struct files*, int, float*);
void terror(char *);
void alert(char *,char *);

void debug_print(char *, int, float*, int);
int  transpose_bin_txt(struct params*, char**);
void qsort_c(float*array,int l,int r,int *index);
int partition( float * a, int l, int r, int *indexes);

// related to Qnorm

void qnorm_b(struct params*, struct files*);


char ** load_probe_id(struct files*, struct params *);

// Load a text-tab-2cols file ans store a bin file
int text_to_bin(struct params *, char**, struct files *);


// ===============================================================================

#endif
