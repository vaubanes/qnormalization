/******************************************************************************
* FILE:qfunc.c
*
* general functions for Qnormalization
*
*  Author : José Manuel Mateos
*
* LAST REVISED: 21/12/09
******************************************************************************/

#include "qfunc.h"
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// This function shows a error and end the program
void terror(char *message) {
  printf(message);
  exit(1);
}

// Read the parameters passed as argument
Params *commandline(int argc, char *argv[]) {
  int i;
  char option;
  Params *parameters;
  parameters=(Params*)malloc(sizeof(Params));

  /* default values-------------------------------------*/
  strcpy(parameters->file_list_experiments,"qInput.txt");
  strcpy(parameters->file_out, "qOut.bin");
  parameters->transpose			 = 0;
  parameters->block_size		 = BLOCK_SIZE;
  parameters->num_genes     	 = DEFAULT_NUM_GENES ;
  parameters->num_experiments    = DEFAULT_NUM_EXPERIMENTS;
  /*----------------------------------------------------------*/

  for (i=1; i<argc; i++)  {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') terror("dash needed at argument (sintx error)");

    option = toupper(argv[i][1]);

    // On-Off flags----------------
    if (option == 'T') {
      parameters->transpose=1;
      continue;
    }

    // -argum=value
    if (argv[i][2]!='=') terror("equal symbol needed at argument (sintx error)");

    switch (option) {
    case 'I':
      strcpy(parameters->file_list_experiments,&argv[i][3]);
      break;
    case 'O':
      strcpy(parameters->file_out, &argv[i][3]);
      break;
    case 'E':
      parameters->num_experiments    = atoi(&argv[i][3]);
      break;
    case 'G':
      parameters->num_genes    = atof(&argv[i][3]);
      break;
    case 'N':
      parameters->block_size   = atoi(&argv[i][3]);
      break;
    default:
      terror("Unknown parameter (see syntax)");
    }
  }
  return parameters;
}



// Load to memory a list of files
// datafile format: fileName[tab]nGenes[tab][format][newLINE]
InfoFile* load_input_files(Params *p) {

  FILE *fileInput;
  InfoFile*info=NULL;
  char line[MAX_SIZE_LINE],line2[MAX_SIZE_LINE],t;
  int index=0,j,g;

  if ((fileInput=fopen(p->file_list_experiments,"rt"))==NULL) terror("opening input file");

  if ((info=(InfoFile*)calloc(p->num_experiments,sizeof(InfoFile)))==NULL)
    terror("memory for list of files");

  fgets(line,MAX_SIZE_LINE,fileInput);
  j=3;
  while (!feof(fileInput) && j==3) {

    if (line[0]!='@') {
      j=sscanf(line,"%s\t%d\t%c\n",line2,&g,&t);
      if (index==p->num_experiments) {
        fprintf(stderr,"[WARNING] more than %d lines... using firts %d as filenames\n",index,index);
        p->num_experiments=index;
        return info;
      }

      if (line2[strlen(line2)-1]=='\t'||line2[strlen(line2)-1]==' ') line2[strlen(line2)-1]=0;

      if (strlen(line2)>0) {
        info[index].file_name=(char*)strdup(line2);
        info[index].num_genes   =g;
        info[index].file_type=t;
        index++;
      }
    }
    fgets(line,MAX_SIZE_LINE,fileInput);

  }
  fclose(fileInput);
  if (index!=p->num_experiments) {
    fprintf(stderr,"[WARNING] only %d files.. nExp=%d\n",index,index);
    p->num_experiments = index;
  }

  return info;
}


// Input returns ordered and Index contains the origial position
int qnorm(double *input, int *dindex, int num_genes) {

  int j;

  for (j=0; j<num_genes;j++) dindex[j]=j; // init the indexes array

  quick_sort(input,0,num_genes-1,dindex); // Quicksort

  return 1;
}


// Add a experiment to average vector
int accumulate_row(Average *average, double *input , int num_genes) {

  int i;

  for (i=0;i<num_genes;i++) {

    average[i].summation+=input[i];
    average[i].num_experiments++;
  }

  return 0;

}


// Calculate number of blocks
int calculate_blocks(int num_experiments, int block_size) {

  int numblocks;
  float numblocks_float = 0;


  if (block_size > num_experiments) {
    numblocks = 1;
  } else {

    numblocks_float = (float) num_experiments / (float) block_size;
  }

  numblocks = (int) ceil(numblocks_float);



  return numblocks;

}

// Calculate inicital blocks
int calculate_initial_blocks(int numblocks, int num_processors) {

  int nblocks;

  if (num_processors  > numblocks) {
    nblocks = numblocks;
  } else {
    nblocks = num_processors;
  }


  return nblocks;
}

// Calculate the index init and end of a block
int calculate_index_blocks(int *index,int num_experiments,int block_size) {

  const int index1 = index[1] + 1;

  int index2 = index1 + block_size -1 ;


  if (index2 > (num_experiments -1)) {
    index2 = num_experiments -1;
  }

  index[0] = index1;
  index[1] = index2;

  return 0;

}


// Sort the values of a experiment
void quick_sort(double *array,int l,int r,int *index) {

  int j;

  if ( l < r ) {
    // divide and conquer
    j = partition( array, l, r,index);
    //  j=(l+r)/2;
    quick_sort( array, l, j-1,index);
    quick_sort( array, j+1, r,index);
  }

}



// Auxiliar function of the quick sort divide the array in two arrays
int partition( double* a, int l, int r, int *indexes) {

  int i=l;
  int j=r+1;
  int k;
  double t;
  double const pivot = a[l];
  //i = l;
  //j = r+1;

  while ( 1) {
    do {
      ++i;
    } while ( a[i] <= pivot && i <= r );
    do {
      --j;
    } while ( a[j] > pivot );
    if ( i >= j ) break;
    t = a[i];
    a[i] = a[j];
    a[j] = t;
    k=indexes[i];
    indexes[i]=indexes[j];
    indexes[j]=k;
  }
  t = a[l];
  a[l] = a[j];
  a[j] = t;
  k=indexes[l];
  indexes[l]=indexes[j];
  indexes[j]=k;
  return j;
}




// Transpose from Disk to Disk the binary matrix into a tab delimited text file
int transpose_matrix(Params *p) {
  FILE *file_out;
  double value, **matrix;
  int i,j;
  char new_name[MAX_SIZE_LINE];
  const int num_genes=p->num_genes;
  const int num_experiments=p->num_experiments;

  if ((file_out=fopen(p->file_out,"rb"))==NULL)
    terror("[Bin2Text] opening binary output file");

  if ((matrix=(double **)calloc(num_genes,sizeof(double*)))==NULL) terror("[Bin2Text] memory for index1");
  for (i=0; i<num_genes;i++)
    if ((matrix[i]=(double *)calloc(num_experiments,sizeof(double)))==NULL) terror("[Bin2Text] memory for index2 full matrix");

  for (i=0; i<num_experiments;i++) {
    for (j=0; j<num_genes; j++) {
      fread(&value, 1, sizeof(double), file_out);
      matrix[j][i]=value;
    }
  }
  fclose(file_out);

  // Save trasposed text tabulated

  sprintf(new_name,"%s.txt",p->file_out);
  if ((file_out=fopen(new_name,"wt"))==NULL)
    terror("[Bin2Text] opening tabulated text file out");

  for (i=0; i<num_genes;i++) {
    fprintf(file_out,"%i\t",i);
    for (j=0; j<num_experiments; j++) {
      if (j) fprintf(file_out,"\t");
      fprintf(file_out,"%lf",matrix[i][j]);
    }
    fprintf(file_out,"\n");
  }

  fclose(file_out);

  return 1;
}



// Load a result of a block calculate for a slave processor
void load_experiment_from_file(InfoFile*flist, int col, double *dataIn, int num_genes ) {

  FILE *file_temp;

  int i;
  double x;

  switch (flist[col].file_type) {
  case 't':
    if ((file_temp=fopen(flist[col].file_name,"rt"))==NULL) terror("opening input file (LoadFIle function)");
    for (i=0;i<num_genes;i++) {
      fscanf(file_temp,"%lf\n",&x);
      dataIn[i]=x;
    }
    fclose(file_temp);
    break;
  default :
    terror("unknow type of file");
  }

}


// DEBUG function for arrays
void debug_print(char *s, double* d, int n) {
  int j;
  fprintf(stderr, "%s------------\n",s);
  for (j=0;j<n;j++) fprintf (stderr,"%f ", d[j]);
  fprintf(stderr,"\n");
}








