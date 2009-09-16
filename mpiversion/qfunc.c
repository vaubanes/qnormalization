/******************************************************************************
* general functions for Qnormalization


* LAST REVISED: 23/02/09
******************************************************************************/

#include "qfunc.h"
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void terror(char *s) {
  printf(s);
  exit(1);
}


Params *commandline(int argc, char *argv[]) {
  int i;
  char c;
  Params *p;
  p=(Params*)malloc(sizeof(Params));

  /* default values-------------------------------------*/
  strcpy(p->flist_experiments,"qInput.txt");
  strcpy(p->fout, "qOut.bin");
  p->traspose		 = 0;
  p->block_size		 = BLOCK_SIZE;
  p->n_genes     	 = NGEN;
  p->n_experiments   = NEXP;
  /*----------------------------------------------------------*/

  for (i=1; i<argc; i++)  {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') terror("dash needed at argument (sintx error)");

    c = toupper(argv[i][1]);

    // On-Off flags----------------
    if (c == 'T') {
      p->traspose=1;
      continue;
    }

    // -argum=value
    if (argv[i][2]!='=') terror("equal symbol needed at argument (sintx error)");

    switch (c) {
    case 'I':
      strcpy(p->flist_experiments,&argv[i][3]);
      break;
    case 'O':
      strcpy(p->fout, &argv[i][3]);
      break;
    case 'E':
      p->n_experiments    = atoi(&argv[i][3]);
      break;
    case 'G':
      p->n_genes    = atof(&argv[i][3]);
      break;
    case 'N':
      p->block_size   = atoi(&argv[i][3]);
      break;
    default:
      terror("Unknown parameter (see syntax)");
    }
  }
  return p;
}



// Load to memory a list of files
// datafile format: fileName[tab]nGenes[tab][format][newLINE]
Files* load_input_files(Params *p) {
  FILE *f;
  Files*L=NULL;
  char line[MAXLIN],lin2[MAXLIN],t;
  int N=0,j,g;

  if ((f=fopen(p->flist_experiments,"rt"))==NULL) terror("opening input file");

  if ((L=(Files*)calloc(p->n_experiments,sizeof(Files)))==NULL)
    terror("memory for list of files");

  fgets(line,MAXLIN,f);
  j=3;
  while (!feof(f) && j==3) {

    if (line[0]!='@') {
      j=sscanf(line,"%s\t%d\t%c\n",lin2,&g,&t);
      if (N==p->n_experiments) {
        fprintf(stderr,"[WARNING] more than %d lines... using firts %d as filenames\n",N,N);
        p->n_experiments=N;
        return L;
      }

      if (lin2[strlen(lin2)-1]=='\t'||lin2[strlen(lin2)-1]==' ') lin2[strlen(lin2)-1]=0;

      if (strlen(lin2)>0) {
        L[N].fname=(char*)strdup(lin2);
        L[N].num_genes   =g;
        L[N].fType=t;
        N++;
      }
    }
    fgets(line,MAXLIN,f);

  }
  fclose(f);
  if (N!=p->n_experiments) {
    fprintf(stderr,"[WARNING] only %d files.. nExp=%d\n",N,N);
    p->n_experiments = N;
  }

  return L;
}


// input returns ordered and Index contains the origial position
int qnorm(double *input, int *dindex, int num_genes) {

  int j;

  for (j=0; j<num_genes;j++) dindex[j]=j; // init the indexes array

  quick_sort(input,0,num_genes-1,dindex); // Quicksort


  return 1;
}

int accumulate_row(Average *average, double *input , int num_genes) {
  int i;

  for (i=0;i<num_genes;i++) {

    average[i].av+=input[i];
    average[i].num++;
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

int calculate_index_blocks(int *index,int num_experiments,int block_size) {



  int index1,index2;

  index1 = index[1] + 1;

  index2 = index1 + block_size -1 ;


  if (index2 > (num_experiments -1)) {
    index2 = num_experiments -1;
  }

  index[0] = index1;
  index[1] = index2;



  return 0;

}




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




int partition( double* a, int l, int r, int *indexes) {
  int i=l;
  int j=r+1;
  int k;
  double t,pivot;
  pivot = a[l];
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


/* DUMMY Functions for initial debuging ==========================*/

// Transpose from Disk to Disk the binary matrix into a tab delimited text file

int transpose_matrix(Params *p) {
  FILE *f;
  double val, **mat;
  int i,j;
  char newname[MAXLIN];
  int num_genes=p->n_genes;
  int num_experiments=p->n_experiments;

  if ((f=fopen(p->fout,"rb"))==NULL)
    terror("[Bin2Text] opening binary output file");

  if ((mat=(double **)calloc(num_genes,sizeof(double*)))==NULL) terror("[Bin2Text] memory for index1");
  for (i=0; i<num_genes;i++)
    if ((mat[i]=(double *)calloc(num_experiments,sizeof(double)))==NULL) terror("[Bin2Text] memory for index2 full matrix");

  for (i=0; i<num_experiments;i++) {
    for (j=0; j<num_genes; j++) {
      fread(&val, 1, sizeof(double), f);
      mat[j][i]=val;
    }
  }
  fclose(f);

  // Save trasposed text tabulated

  sprintf(newname,"%s.txt",p->fout);
  if ((f=fopen(newname,"wt"))==NULL)
    terror("[Bin2Text] opening tabulated text file out");

  for (i=0; i<num_genes;i++) {
    fprintf(f,"%i\t",i);
    for (j=0; j<num_experiments; j++) {
      if (j) fprintf(f,"\t");
      fprintf(f,"%lf",mat[i][j]);
    }
    fprintf(f,"\n");
  }

  fclose(f);

  return 1;
}


// =======================================================================
// Dummy functions : read a datafile (case t: text: one line/one value)

void load_parcial_result(Files*flist, int col, double *dataIn, int num_genes ) {

  FILE *f;

  int i;
  double x;


  switch (flist[col].fType) {
  case 't':
    if ((f=fopen(flist[col].fname,"rt"))==NULL) terror("opening input file (LoadFIle function)");
    for (i=0;i<num_genes;i++) {
      fscanf(f,"%lf\n",&x);
      dataIn[i]=x;
    }
    fclose(f);
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








