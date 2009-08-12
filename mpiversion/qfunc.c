/******************************************************************************
* general functions for Qnormalization


* LAST REVISED: 23/02/09
******************************************************************************/

#include "qfunc.h"
#include <string.h>
#include <ctype.h>


void terror(char *s) {
  printf(s);
  exit(1);
}


struct params *commandline(int argc, char *argv[]) {
  int i;
  char c;
  struct params *p;
  p=(struct params*)malloc(sizeof(struct params));

  /* default values-------------------------------------*/
  strcpy(p->fListName,"qInput.txt");
  strcpy(p->fOutName, "qOut.bin");
  p->Traspose= 0;
  p->MemIndex= 1;
  p->nP      = MAXnP;
  p->nG      = NGEN;
  p->nE      = NEXP;
  /*----------------------------------------------------------*/

  for (i=1; i<argc; i++)  {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') terror("dash needed at argument (sintx error)");

    c = toupper(argv[i][1]);

    // On-Off flags----------------
    if (c == 'T') {
      p->Traspose=1;
      continue;
    }
    if (c == 'D') {
      p->MemIndex=0;
      continue;
    }

    // -argum=value
    if (argv[i][2]!='=') terror("equal symbol needed at argument (sintx error)");

    switch (c) {
    case 'I':
      strcpy(p->fListName,&argv[i][3]);
      break;
    case 'O':
      strcpy(p->fOutName, &argv[i][3]);
      break;
    case 'E':
      p->nE    = atoi(&argv[i][3]);
      break;
    case 'G':
      p->nG    = atof(&argv[i][3]);
      break;
    case 'N':
      p->nP   = atoi(&argv[i][3]);
      break;
    default:
      terror("Unknown parameter (see syntax)");
    }
  }
  return p;
}



// Load to memory a list of files
// datafile format: fileName[tab]nGenes[tab][format][newLINE]
struct Files* load_input_files(struct params *p) {
  FILE *f;
  struct Files*L=NULL;
  char line[MAXLIN],lin2[MAXLIN],t;
  int N=0,j,g;

  if ((f=fopen(p->fListName,"rt"))==NULL) terror("opening input file");

  if ((L=(struct Files*)calloc(p->nE,sizeof(struct Files)))==NULL)
    terror("memory for list of files");

  fgets(line,MAXLIN,f);
  j=3;
  while (!feof(f) && j==3) {

    if (line[0]!='@') {
      j=sscanf(line,"%s\t%d\t%c\n",lin2,&g,&t);
      if (N==p->nE) {
        fprintf(stderr,"[WARNING] more than %d lines... using firts %d as filenames\n",N,N);
        p->nE=N;
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
  if (N!=p->nE) {
    fprintf(stderr,"[WARNING] only %d files.. nExp=%d\n",N,N);
    p->nE = N;
  }

  return L;
}


// input returns ordered and Index contains the origial position
int qnorm(double *input, int *dindex, int num_genes) {

  int j;

  for (j=0; j<num_genes;j++) dindex[j]=j; // init the indexes array

  sort(input,0,num_genes-1,dindex); // Quicksort


  return 1;
}

int accumulate_row(struct Average *average, double *input , int num_genes) {
  int i;

  for (i=0;i<num_genes;i++) {

	average[i].Av+=input[i];
	average[i].num++;
  }

  return 0;

}



// Calculate number of blocks
int calculate_blocks(int num_experiments) {

  int numblocks;
  float numblocks_float = 0;


  if (BLOCK_SIZE > num_experiments) {
    numblocks = 1;
  } else {

    numblocks_float = (float) num_experiments / (float) BLOCK_SIZE;
  }

  //numBlocks = (int) ceil(numBlocksF);

  numblocks = (int) numblocks_float;

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

int calculate_index_blocks(int *index,int num_experiments) {

  int index1,index2;

  index1 = index[1] + 1;

  index2 = index1 + BLOCK_SIZE -1 ;

  if (index2 > (num_experiments -1)) {
    index2 = num_experiments -1;
  }

  index[0] = index1;
  index[1] = index2;

  return 0;

}




void sort(double *array,int l,int r,int *index) {
  int j;
  if ( l < r ) {
    // divide and conquer
    j = partition( array, l, r,index);
    //  j=(l+r)/2;
    sort( array, l, j-1,index);
    sort( array, j+1, r,index);
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

int transpose_matrix(struct params *p) {
  FILE *f;
  double val, **mat;
  int i,j;
  char NewName[MAXLIN];
  int nG=p->nG;
  int nE=p->nE;

  if ((f=fopen(p->fOutName,"rb"))==NULL)
    terror("[Bin2Text] opening binary output file");

  if ((mat=(double **)calloc(nG,sizeof(double*)))==NULL) terror("[Bin2Text] memory for index1");
  for (i=0; i<nG;i++)
    if ((mat[i]=(double *)calloc(nE,sizeof(double)))==NULL) terror("[Bin2Text] memory for index2 full matrix");

  for (i=0; i<nE;i++) {
    for (j=0; j<nG; j++) {
      fread(&val, 1, sizeof(double), f);
      mat[j][i]=val;
    }
  }
  fclose(f);

  // Save trasposed text tabulated

  sprintf(NewName,"%s.txt",p->fOutName);
  if ((f=fopen(NewName,"wt"))==NULL)
    terror("[Bin2Text] opening tabulated text file out");

  for (i=0; i<nG;i++) {
    fprintf(f,"%i\t",i);
    for (j=0; j<nE; j++) {
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

void load_parcial_result(struct Files*flist, int col, double *dataIn, int num_genes ) {

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








