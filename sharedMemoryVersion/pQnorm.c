/******************************************************************************
* FILE: pQnorm.c
* DESCRIPTION:
*   parallel version
*   Improved sequential prototype for Qnorm 
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: O.Trelles (23 Feb.09)
* 23.Feb.09  : Using command line argums
*              Qnorm [-o=Value]  see below Command line params
*
*              fList  : contains a list of nExp filenames (gene-expression data files) 
*                          line format: fileName[TAB]nGenes[TAB]FileType[NEWLINE]
*              nRows  : number of genes in each file 
*              Normalised.fname where the normalised values will be stored 
*                       (as a text-tabulated matrix)
*              mode : m: keep the Index matrix in memory ( d: in disk)
*                       THIS VERSION REQUIRES ALL FILES WITH THE SAME NUMBER OF GENES
*
*
*  Command line Parameters
*         sintaxis:    Qnorm [-Option=value]... 
*
*  Option  Description              Default value     alternative Values
*  ------  -----------------------  -------------     ------------------
*  -p      Number of Processors        4               Only in parallel version
*  -i      File name (list of files)   qInput.txt      valid existing pathname
*  -o      Output binary matrix        qOut.bin        binary by columns file
*  -e      Number of experiments       2               positive integer
*  -g      NUmber of genes             15              positive integer
*  -M      Index Matrix in mem         D (in disk)     -M (in memory)
*  -V      Verbose mode                Not             -V  
* ---------------------------------------------------------------------------
*  -t      Traspose the fileOut        Not             -T (yes)
   NOT AVAILABLE FROM version 5
*
*	@returns >0 if everything was fine <0 if there was an error

   pQnorm3  : since the memory allocation is performed only by thread 0
              the starting point of the parallel section has been moved
              to the init of the main process----
   pQnorm5 : to share a file each thread uses their own handle
             and the file is opened by each thread

             new GeneExpression file format
             probeID[tab]ExpValueIntensity

   pQnorm6 : using quicksort interative (instead of recursive)
             to avoid stack problems (more than 6M points)
*
* LAST REVISED: 27/02/09
******************************************************************************/
#include <omp.h>
#include "Qfunc.c"
#include "quicksort/quicksort.c"


int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;
        p = CommandLine(ac,av);

	if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");
        QNormMain(p,fList);	

	
	return 1;
}


void QNormMain(struct params *p, struct Files* fList){
   double **dataIn;
   int **dIndex;
   double *dataOut;
   int **mIndex;
   int *dInde2;
   struct Average **AvG; // global Average by row (now a matrx)
   int i,ii,j,k;
   int nG=p->nG;
   int nE=p->nE;
   int nP, nP1, tid;
   int From, To, Range; 
   double checkVal=0;

// the partial AvG array used by each thread to get the 
// accumulated value, will be shared using a matrix (nPxAvGsize)
// to facilitate final global average

    tid = omp_get_thread_num();
    nP = p->nP;  // omp_get_num_threads();
    nP1=nP+1;
    
    if (p->Verbose) fprintf(stderr,"Number of threads = %d\n", nP);

    // Memory===========================================
    // Index array 
    if (p->MemIndex) { // in memory - full
         if ((mIndex=(int **)calloc(nG,sizeof(int*)))==NULL) 
              terror("memory for index1");
         for (ii=0; ii<nG;ii++)
           if((mIndex[ii]=(int *)calloc(nE,sizeof(int)))==NULL) 
             terror("memory for index2 full matrix");
    }
    if ((AvG =(struct Average **)calloc(nP1,sizeof(struct Average*)))==NULL)
       terror("memory for Average array");
    for (ii=0; ii<nP1;ii++)
       if((AvG[ii]=(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) 
              terror("memory for average 2 full matrix");

   // This will always be necessary to decuple the function
    if ((dIndex =(int**)calloc(nP,sizeof(int *)))==NULL)
       terror("memory for dIndex mat");
    for (ii=0; ii<nP;ii++)
       if((dIndex[ii]=(int *)calloc(nG,sizeof(int)))==NULL) 
           terror("memory for dIndex");
    if ((dataIn =(double**)calloc(nP,sizeof(double *)))==NULL)
       terror("memory for dataIn mat");
    for (ii=0; ii<nP;ii++)
       if ((dataIn[ii]=(double *)calloc(nG,sizeof(double)))==NULL) 
          terror("memory for dataIn array");

   // oputput file (by cols) (to share the handle)------
#pragma omp parallel shared(nG, nE,mIndex,dataIn, dIndex, AvG, fList, p) private(i,j,k,tid,nP1, dataOut, From, To, Range)
 { // Open General parallel section [0]

    FILE *fI;
    long posbase;
    int  posbasei;
    int pos;


    tid = omp_get_thread_num();
    nP = omp_get_num_threads();
    if (nP != p->nP) terror("something wrong in nP");

   if (!p->MemIndex) { // master opens the file
      if ((fI=fopen("~tmp","wb"))==NULL) terror("opening tmp-index file");
   }
   // LOAD DISTRIBUTION------------------
   Range = nE / nP;
   From = tid * Range;
   To   = (tid+1)*Range;
   if (To > nE) To = nE;

   // each P initialize
   for (j=0; j< nG;j++) { // init Accumulation array 
      AvG[tid][j].Av=0;        // =HUGE_VAL; ???
      AvG[tid][j].num=0;
      if (tid==0) {
        AvG[nP][j].Av=0;        // =HUGE_VAL; ???
        AvG[nP][j].num=0;
      }
   }

   // QNORM ===============================================================
   if (p->Verbose) fprintf(stderr,"[1st-Step]");

   for (i=From; i< To; i++) { // Qnorm for each datafile: STEP 1
        LoadFile(fList, i, dataIn[tid]);
        if (p->Verbose) { fprintf(stderr,"."); fflush(stderr);}

#ifdef DEBUG
        DebugPrint("Load",tid, dataIn[tid], fList[i].nG); 
#endif
        // dataIn returns ordered and Index contains the origial position
        Qnorm1(dataIn[tid], dIndex[tid], fList[i].nG); 

#ifdef DEBUG
        DebugPrint("Sorted",tid, dataIn[tid], nG);         
#endif
        
        AccumulateRow(AvG[tid], dataIn[tid] , nG);

        // now decide how to proceed with indexes
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) 
            mIndex[j][i]= dIndex[tid][j];
        } else {         // in disk
          pos=fList[i].pos;
          fseek(fI, nG*pos*sizeof(int), SEEK_SET);
          fwrite(dIndex[tid], sizeof(int), nG, fI);
        }

#ifdef DEBUG
        fprintf(stderr,"Index tid=%d (col=%d)\n",tid,i);
        for (j=0;j<nG;j++) fprintf (stderr,"%d ", dIndex[tid][j]); fprintf(stderr,"\n");
#endif


   } // end "for i" Qnorm for each datafile: STEP 1

   if (!p->MemIndex)  fclose(fI);
}

   // HERE only one thread
   // Row average     [use col=tid=0] ------------------------------

   if (tid==0) {
     for (i=0;i<nG;i++) {
       for (j=0;j<nP;j++){
         AvG[nP][i].Av +=AvG[j][i].Av;
         AvG[nP][i].num+=AvG[j][i].num;
       }
       AvG[nP][i].Av /=AvG[nP][i].num;
       checkVal +=AvG[nP][i].Av;
      }
      checkVal /=nG;
      if (p->Verbose) fprintf(stderr, "checkVal=%lf\n",checkVal);
      fprintf(stderr, "tid=%d checkVal=%lf\n",tid,checkVal); fflush(stderr);

      // Now copy the global averge into the local ones
      // *distribute* the global array
      for (i=0;i<nG;i++) {
        for (j=0;j<nP;j++) {
           AvG[j][i].Av =AvG[nP][i].Av;
           AvG[j][i].num=AvG[nP][i].num;
        }
      }


#ifdef DEBUG
     fprintf(stderr, "Row Average----[tid==%d]--------\n",tid);
     for (j=0;j<nG;j++) fprintf (stderr,"%f (%d) ", AvG[tid][j].Av,AvG[tid][j].num); 
     fprintf(stderr,"\n");
#endif

    } // end tid=0



// second step----------------------------------------------------
// sychro needed**

#pragma omp parallel shared(nG, nE,mIndex,dataIn, dIndex, AvG, fList, p) private(i,j,k,tid,nP1, dataOut, From, To, Range)
 { // Open General parallel section [0]

    FILE *fI2, *fOut;
    long posbase;
    int  posbasei;
    int pos;

    tid = omp_get_thread_num();
    nP = omp_get_num_threads();
    if (nP != p->nP) terror("something wrong in nP");

   // LOAD DISTRIBUTION------------------
   Range = nE / nP;
   From = tid * Range;
   To   = (tid+1)*Range;
   if (To > nE) To = nE;

   if ((fOut=fopen(p->fOutName,"wb"))==NULL) {
      fprintf(stderr,"file : %s\n",p->fOutName); fflush(stderr);
       terror("opening OUTPUT file");
   }
   // Finally produce the ORDERED output file [STEP 2]--------------------

   if (!p->MemIndex) { 
      if ((fI2=fopen("~tmp","rb"))==NULL) 
         terror("opening tmp-index for reading file");
   }

   if ((dataOut=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataOut array");
   if (p->Verbose) { fprintf(stderr,"\n[2nd Step]"); fflush(stderr);}
   if((dInde2=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

   if (p->Verbose) {
      fprintf(stderr, "[tid=%d] 2nd Step\n",tid); fflush(stderr);
   }
   posbase = (long)nG*sizeof(double);
   posbasei= (long)(nG)*sizeof(int);
   for (i=From;i<To;i++) {
        if (p->MemIndex) { // in memory - full
          for (j=0;j<nG;j++) { 
            dInde2[j]=mIndex[j][i];
          }
        } else {
          pos = fList[i].pos;
          // fprintf(stderr, "[tid=%d] lee pos=%d(i=#file=%d) posbasei=%d\n",tid,pos,i,posbasei*pos); fflush(stderr);
          fseek(fI2, posbasei*pos, SEEK_SET);
          fread(dInde2, sizeof(int), nG, fI2);
        }

        // complete the output vector
        for (j=0;j<nG;j++) 
          dataOut[dInde2[j]]=AvG[tid][j].Av; // OJO

        pos = fList[i].pos;
        fseek(fOut, posbase*pos, SEEK_SET);
        fwrite(dataOut, sizeof(double), nG, fOut);
        //fprintf(stderr, "[tid=%d] ESC pos=%ld\n",tid,posbase*pos); fflush(stderr);
           
   }
   fclose(fOut);
   if (!p->MemIndex) fclose(fI2);

} // close General parallel section [0]

    if (p->Verbose) fprintf(stderr,"done...\n");   

    return ;

} // end MAIN------------------


// input returns ordered and Index contains the origial position

int Qnorm1(double *input, int *dIndex, int nG){
	int i,j,k,n;

	for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array


        quicksort(input,dIndex,nG); // iterative quicksort (stack problems using recursivity)
//	QsortC(input,0,nG-1,dIndex); // Quicksort 
        return 1;
}

void AccumulateRow(struct Average *AvG, double *input , int nG){
        int i;
        double tot=0;
 	
        for (i=0;i<nG;i++) {

		AvG[i].Av+=input[i];
		AvG[i].num++;
                tot+=input[i];
	}

	return;

}

