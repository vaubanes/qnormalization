/******************************************************************************
* general functions for Qnormalization


* LAST REVISED: 23/02/09
******************************************************************************/

#include "Qfunc.h"



struct params *CommandLine(int argc, char *argv[])
{
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
  p->Verbose = 0;
  p->GPU = -1;
  p->threads = 256;
  p->blocks = 256;
  p->shared = 1024;

  /*----------------------------------------------------------*/

  for (i=1; i<argc; i++)  {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') terror("dash needed at argument (sintx error)");

    c = toupper(argv[i][1]);

    // On-Off flags----------------
//    if (c == 'T') {p->Traspose=1; continue;}
//    if (c == 'D') {p->MemIndex=0; continue;}
    if (c == 'V') {p->Verbose =1; continue;}

    // -argum=value
    if (argv[i][2]!='=') terror("equal symbol needed at argument (sintx error)");

    switch(c)
    {
      case 'I': strcpy(p->fListName,&argv[i][3]);  break;
      case 'O': strcpy(p->fOutName, &argv[i][3]);  break;
      case 'E': p->nE    = atoi(&argv[i][3]);      break;
      case 'G': p->nG    = atof(&argv[i][3]);      break;
      case 'M': 1;                                 break;
      case 'P': p->nP   = atoi(&argv[i][3]);       break;
      case 'H': p->GPU   = atoi(&argv[i][3]);       break;
      case 'S': p->shared   = atoi(&argv[i][3]);       break;
      case 'B': p->blocks   = atoi(&argv[i][3]);       break;
      case 'T': p->threads   = atoi(&argv[i][3]);       break;
      default: terror("Unknown parameter (see syntax)");
    }
  }
  return p;
}


// panic message-----------------
void terror(char *s)
{ fprintf(stdout,"\n[ERR**] %s **\n",s); exit(-1); }

// warning message---------------
void Alerta(char *s,char *s1)
{ fprintf(stdout,"\n[WARNING] %s : %s\n",s,s1); }
  

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
        while(!feof(f)) {

           if (line[0]!='@') {

             j=sscanf(line,"%s\t%d\t%c\n",lin2,&g,&t);         
             if (N==p->nE) {
                if (p->Verbose) 
                   fprintf(stderr,"[WARNING] more than %d lines in file LIST... using first %d as Experiments\n",N,N);
                p->nE=N;
                return L;
             }
         
             if (lin2[strlen(lin2)-1]=='\t'||lin2[strlen(lin2)-1]==' ') lin2[strlen(lin2)-1]=0;

             if (strlen(lin2)>0) {
               L[N].fname=(char*)strdup(lin2);
               L[N].nG   =g;
               L[N].fType=t;
               L[N].pos=N;
               if (g > p->nG) {
                  if (p->Verbose)
                    Alerta("more genes in file than in parameter, work with min value",L[N].fname);
                  L[N].nG   =p->nG;
               }
               if (g < p->nG) {
                  if (p->Verbose)
                    Alerta("more genes in parametern than in file, work with min value",L[N].fname);
                  p->nG = L[N].nG;
               }

               N++;
             }
           }
           fgets(line,MAXLIN,f);

       }
       fclose(f);
       if (N!=p->nE) {
           if (p->Verbose) fprintf(stderr,"[WARNING] only %d files.. nExp=%d\n",N,N);
           p->nE = N;
         }

       return L;
}

 

void QsortC(double *array,int l,int r,int *index)
{
   int j;
   if( l < r ) {
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
		
   while( 1) {
     do{ ++i; 
      }while( a[i] <= pivot && i <= r );
     do{ --j; 
      }while( a[j] > pivot );
     if( i >= j ) break;
     t = a[i]; a[i] = a[j]; a[j] = t;
     k=indexes[i];indexes[i]=indexes[j];indexes[j]=k;
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   k=indexes[l];indexes[l]=indexes[j];indexes[j]=k;
   return j;
}


/* DUMMY Functions for initial debuging ==========================*/

// Transpose from Disk to Disk the binary matrix into a tab delimited text file

int TransposeBin2Txt(struct params *p,char**probeID){
	FILE *f,*f1;
        float **mat,val;
	int i,j,k;
        char NewName[MAXLIN];
        int nG=p->nG;
        int nE=p->nE;
        int BLQ=1000;
        int From=0,pos;
    
        if (BLQ>nG) BLQ=nG;

        if ((f=fopen(p->fOutName,"rb"))==NULL)
           terror("[Bin2Text] opening binary output file");

        sprintf(NewName,"%s.txt",p->fOutName);
        if ((f1=fopen(NewName,"wt"))==NULL) 
           terror("[Bin2Text] opening tabulated text file out");

        if ((mat=(float **)calloc(nE,sizeof(float*)))==NULL)
            terror("[Bin2Text] memory for index1");
        for (i=0; i<nE;i++)
          if((mat[i]=(float *)calloc(BLQ,sizeof(float)))==NULL) 
              terror("[Bin2Text] memory for index2 full matrix");
  
        while (From<=nG) { 
	  for (i=0; i<nE;i++){ // read a BLQ of genes from each Exp
             pos=(i*nG+From)*sizeof(float);
             fseek(f, pos, SEEK_SET); 
             fread(mat[i], BLQ, sizeof(float), f);
/*---
             for(k=0;k<BLQ;k++) { 
               fread(&val, 1, sizeof(double), f);
               mat[i][k]=val;
               //printf("lee=%lf ",val);
             }
---*/
          }
          // Save a BLQ trasposed text tabulated

	  for (i=0; i<BLQ&&From+i<nG;i++){
             fprintf(f1,"%s",probeID[From+i]);
	     for (j=0; j<nE; j++)
                  fprintf(f1,"\t%9.2f",mat[j][i]);
             fprintf(f1,"\n");
           }
           From+=BLQ;
        }
	
        fclose(f);
        fclose(f1);
	
	return 1;
}


// =======================================================================
// Dummy functions : read a datafile (case t: text: one line/one value)

void LoadFile(struct Files*fList, int col, unsigned int *dataIn){
        FILE *f;
        char line[MAXLIN],lin2[MAXLIN],t;
        int N=0,j,g,i;
        float x;
        
        switch(fList[col].fType){
           case 't':
                     if ((f=fopen(fList[col].fname,"rt"))==NULL) {
                        fprintf(stderr,"file:%s\n",fList[col].fname);
                        terror("opening input file (LoadFIle function)");
                     }
                     for (i=0;i<fList[col].nG;i++) {
                       fscanf(f," %*s %f ",&x);
                       dataIn[i]=(unsigned int) x;
                     }
                     fclose(f);
                     break;
           default : terror("unknow (LoadFile) type of file");
        }

}

char ** LoadprobeID(struct Files* fList, struct params *p){
        FILE *f;
        int i;
        char probeID[255],**probe;
        int col=0;
        float x;

        if ((probe=(char **)calloc(p->nG,sizeof(char*)))==NULL) 
           terror("[Bin2Text] memory for probes");
        switch(fList[col].fType){
           case 't':
                     if ((f=fopen(fList[col].fname,"rt"))==NULL) {
                        fprintf(stderr,"file:%s\n",fList[col].fname);
                        terror("opening input file (LoadFIle function)");
                     }
                     for (i=0;i<fList[col].nG;i++) {
                       fscanf(f,"%s\t%f\n",probeID,&x);
                       probe[i]=strdup(probeID);
                     }
                     fclose(f);
                     break;
           default : terror("unknow (loadProbeID) type of file");
        }
       return probe;

}


// DEBUG function for arrays
void DebugPrint(char *s, int tid,float* d, int n){
     int j;
     fprintf(stderr, "%s----[tid=%d]--------\n",s,tid);
     for (j=0;j<n;j++) fprintf (stderr,"%4.1f ", d[j]); 
     fprintf(stderr,"\n");
}
