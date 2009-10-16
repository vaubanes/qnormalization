/******************************************************************************
* FILE: pQnorm.c
* DESCRIPTION:
*   parallel version
*   Improved sequential prototype for Qnorm 
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*  
*   pQnorm7  (a) read InFiles, sort and store in disk the index matrix
*               Produces the "Average Vector" file Out
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
*  NOT AVAILABLE FROM version 5
*
*  -t      Traspose the fileOut        Not             -T (yes)
*          Deprecated in v6. Transposing procedure is performed by Qtrans.c
*  -M      Index Matrix in mem         D (in disk)     -M (in memory)
*          Deprecated in v7. Using only Disk version
*          (code for managing the Index matrix in memory is still
*           valid to keep compatibility con previous versions)
*
*	@returns >0 if everything was fine <0 if there was an error
*
*  pQnorm3  : since the memory allocation is performed only by thread 0
*             the starting point of the parallel section has been moved
*             to the init of the main process----
*  pQnorm5 : to share a file each thread uses their own handle
*            and the file is opened by each thread
*            new GeneExpression file format
*            probeID[tab]ExpValueIntensity
*  pQnorm6 : using quicksort interative (instead of recursive)
*            to avoid stack problems (more than 6M points)
*  pQnorm7  (a) read InFiles, sort and store in disk the index matrix
*               Produces the "Average Vector" file Out
*           (b) read index and AvVector and produce the Qnorm matrix stored by cols
*
* LAST REVISED: 13/05/09
******************************************************************************/
#include <omp.h>
#include "qnorm_a.h"


int main(int ac, char **av){

    struct Files *file_list=NULL;
    struct params *parameters=NULL;

    parameters = command_line(ac,av);

	if ((file_list=load_list_files(parameters))==NULL)
           terror("Loading list of files");

	parameters->MemIndex=0; // working in disk
    qnorm_main(parameters,file_list);
	
	return 1;
}


void qnorm_main(struct params *parameters, struct Files* fList){

   float**dataIn;
   int **dIndex;
   int **mIndex;
   struct Average **average; // global Average by row (now a matrx)
   int i,ii,j;
   int num_genes=parameters->nG;
   int num_experiments=parameters->nE;
   int numProcesors, nP1, tid;
   int from, to, range; 
   float checkVal=0;
   // the partial AvG array used by each thread to get the
   // accumulated value, will be shared using a matrix (nPxAvGsize)
   // to facilitate final global average

    tid = omp_get_thread_num();
    numProcesors = parameters->nP;  // omp_get_num_threads();
    nP1=numProcesors+1;
    
    if (parameters->Verbose) fprintf(stderr,"Number of threads = %d\n", numProcesors);

    // Memory===========================================
    // Index array 
    if (parameters->MemIndex) { // in memory - full
         if ((mIndex=(int **)calloc(num_genes,sizeof(int*)))==NULL)
              terror("memory for index1");
         for (ii=0; ii<num_genes;ii++)
           if((mIndex[ii]=(int *)calloc(num_experiments,sizeof(int)))==NULL)
             terror("memory for index2 full matrix");
    }
    if ((average =(struct Average **)calloc(nP1,sizeof(struct Average*)))==NULL)
       terror("memory for Average array");
    for (ii=0; ii<nP1;ii++)
       if((average[ii]=(struct Average *)calloc(num_genes,sizeof(struct Average)))==NULL)
              terror("memory for average 2 full matrix");

   // This will always be necessary to decuple the function
    if ((dIndex =(int**)calloc(numProcesors,sizeof(int *)))==NULL)
       terror("memory for dIndex mat");
    for (ii=0; ii<numProcesors;ii++)
       if((dIndex[ii]=(int *)calloc(num_genes,sizeof(int)))==NULL)
           terror("memory for dIndex");
    if ((dataIn =(float**)calloc(numProcesors,sizeof(float*)))==NULL)
       terror("memory for dataIn mat");
    for (ii=0; ii<numProcesors;ii++)
       if ((dataIn[ii]=(float*)calloc(num_genes,sizeof(float)))==NULL)
          terror("memory for dataIn array");

   // oputput file (by cols) (to share the handle)------
#pragma omp parallel shared(num_genes, num_experiments,mIndex,dataIn, dIndex, average, fList, parameters) private(i,j,tid,nP1, from, to, range)
 { // Open General parallel section [0]

    FILE *fI;
    int pos;


    tid = omp_get_thread_num();
    numProcesors = omp_get_num_threads();
    if (numProcesors != parameters->nP) terror("something wrong in nP");

   if (!parameters->MemIndex) { // master opens the file
      if ((fI=fopen("~tmp","wb"))==NULL) terror("opening tmp-index file");
   }
   // LOAD DISTRIBUTION------------------
   range = num_experiments / numProcesors;
   from = tid * range;
   to   = (tid+1)*range;
   if (to > num_experiments) to = num_experiments;

   // each P initialize
   for (j=0; j< num_genes;j++) { // init Accumulation array
      average[tid][j].Av=0;   
      average[tid][j].num=0;
      if (tid==0) {
        average[numProcesors][j].Av=0; 
        average[numProcesors][j].num=0;
      }
   }

   // QNORM ===============================================================
   if (parameters->Verbose) fprintf(stderr,"[1st-Step]");

   for (i=from; i< to; i++) { // Qnorm for each datafile: STEP 1
        LoadFile(fList, i, dataIn[tid]);
        if (parameters->Verbose) { fprintf(stderr,"."); fflush(stderr);}

#ifdef DEBUG
        DebugPrint("Load",tid, dataIn[tid], fList[i].nG); //debug example
#endif
        // dataIn returns ordered and Index contains the origial position
        qnorm_sort(dataIn[tid], dIndex[tid], fList[i].nG);

        accumulate_row(average[tid], dataIn[tid] , num_genes);

        // now decide how to proceed with indexes
        if (parameters->MemIndex) { // in memory - full
          for (j=0;j<num_genes;j++)
            mIndex[j][i]= dIndex[tid][j];
        } else {         // in disk :active Option in v7
          pos=fList[i].pos;
          fseek(fI, num_genes*pos*sizeof(int), SEEK_SET);
          fwrite(dIndex[tid], sizeof(int), num_genes, fI);
        }


   } // end "for i" Qnorm for each datafile: STEP 1

   if (!parameters->MemIndex)  fclose(fI);
}

   // HERE only one thread
   // Row average     [use col=tid=0] ------------------------------

   fprintf(stderr, "tid=%d AvG(0)(0)=%f num=%d\n",tid,average[0][0].Av, average[0][0].num);
   if (tid==0) {
     for (i=0;i<num_genes;i++) {
       for (j=0;j<numProcesors;j++){
         average[numProcesors][i].Av +=average[j][i].Av;
         average[numProcesors][i].num+=average[j][i].num;
       }
       average[numProcesors][i].Av /=average[numProcesors][i].num;
       checkVal +=average[numProcesors][i].Av;
      }
      checkVal /=(float)num_genes;
      if (parameters->Verbose) fprintf(stderr, "checkVal=%f\n",checkVal);
      fprintf(stderr, "tid=%d checkVal=%f\n",tid,checkVal); fflush(stderr);

      { // New section v7--------------------
 
        char nfname[500];
        FILE *ff;
        sprintf(nfname, "%s.avg",parameters->fOutName);

        if ((ff=fopen(nfname,"wb"))==NULL) {
           fprintf(stderr,"file : %s\n",nfname); 
           terror("opening OUTPUT file for Average Vector");
        }

        // Now copy the global averge into the local ones
        // *distribute* the global array
          for (j=0;j<num_genes;j++)
             dataIn[0][j] =average[numProcesors][j].Av;
          fwrite(dataIn[0], sizeof(float), num_genes, ff);
          fclose(ff);
     } // end new section (store Average Vector)
   } // end tid=0

    if (parameters->Verbose) fprintf(stderr,"done...\n");   

    return ;

} // end MAIN------------------


// input returns ordered and Index contains the origial position
int qnorm_sort(float*input, int *dIndex, int nG){
	int j;

	for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array

    quicksort(input,dIndex,nG); // iterative quicksort (stack problems using recursivity)

    return 1;
}

void accumulate_row(struct Average *AvG, float*input , int nG){
        int i;
        float tot=0;
 	
        for (i=0;i<nG;i++) {

		AvG[i].Av+=input[i];
		AvG[i].num++;
                tot+=input[i];
	}

	return;

}

