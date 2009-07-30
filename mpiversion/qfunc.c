/******************************************************************************
* general functions for Qnormalization


* LAST REVISED: 23/02/09
******************************************************************************/

#include "qfunc.h"


void terror(char *s)
{
  printf(s);
  exit(1);
}


struct params *CommandLine(int argc, char *argv[]) {
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

/**
// panic message-----------------
void terror(char *s)
{ fprintf(stdout,"\n[ERR**] %s **\n",s); exit(-1); }
**/
// warning message---------------
void Alerta(char *s,char *s1) {
  fprintf(stdout,"\n[WARNING] %s : %s\n",s,s1);
}


// Load to memory a list of files
// datafile format: fileName[tab]nGenes[tab][format][newLINE]

struct Files* LoadListOfFiles(struct params *p) {
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
        L[N].nG   =g;
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
int Qnorm1(double *input, int *dIndex, int nG) {

  int j;

  for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array

  QsortC(input,0,nG-1,dIndex); // Quicksort
  //quicksort(input,dIndex,nG);

  return 1;
}

int AccumulateRow(struct Average *AvG, double *input , int nG) {
  int i;

  for (i=0;i<nG;i++) {

    AvG[i].Av+=input[i];
    AvG[i].num++;
  }

  return 0;

}



// Calculate number of blocks
int calculateBlocks(int nE) {

  int numBlocks;
  float numBlocksF = 0;


  if (BLOCK_SIZE > nE) {
    numBlocks = 1;
  } else {

    numBlocksF = (float) nE / (float) BLOCK_SIZE;
  }

  numBlocks = (int) ceil(numBlocksF);

  return numBlocks;

}

// Calculate inicital blocks
int calculateInitialBlocks(int numBlocks, int nProcesors) {

  int nBlocks;

  if (nProcesors  > numBlocks) {
    nBlocks = numBlocks;
  } else {
    nBlocks = nProcesors;
  }


  return nBlocks;
}

int calculateIndexBlocks(int *index,int nE) {

  int index1,index2;

  index1 = index[1] + 1;

  index2 = index1 + BLOCK_SIZE -1 ;

  if (index2 > (nE -1)) {
    index2 = nE -1;
  }

  index[0] = index1;
  index[1] = index2;

  return 0;

}




void QsortC(double *array,int l,int r,int *index) {
  int j;
  if ( l < r ) {
    // divide and conquer
    j = partition( array, l, r,index);
    //  j=(l+r)/2;
    QsortC( array, l, j-1,index);
    QsortC( array, j+1, r,index);
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

int TransposeBin2Txt(struct params *p) {
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

void LoadFile(struct Files*fList, int col, double *dataIn, int nG ) {

  FILE *f;
  char line[MAXLIN],lin2[MAXLIN],t;
  int N=0,j,g,i;
  double x;


  switch (fList[col].fType) {
  case 't':
    if ((f=fopen(fList[col].fname,"rt"))==NULL) terror("opening input file (LoadFIle function)");
    for (i=0;i<nG;i++) {
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
void DebugPrint(char *s, double* d, int n) {
  int j;
  fprintf(stderr, "%s------------\n",s);
  for (j=0;j<n;j++) fprintf (stderr,"%f ", d[j]);
  fprintf(stderr,"\n");
}








