/******************************************************************************
* FILE: pQnorm.c (v7b)
*
* Important: 
* from v7 the pQnorm procedure was splited in two parts
* part (a) load/sort/storeIndex and finally store the Average Vector (PARALLEL)
* part (b) read index / prepare Qnorm vector / Store (SEQUENTIAL)
* (remember that matrix transposing is performed by Qtrans)
*
*

* DESCRIPTION:
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
*
* LAST REVISED: 13/05/09
******************************************************************************/
#include "Qfunc.c"

int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;
        p = CommandLine(ac,av);
	if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");
        Q2(p,fList);	

	
	return 1;
}

/*-------------------------------------------
This section can be runned in sequential mode (one PE)
since it is a I/O bounded procedure.
However, I will maintain the parallel functionality to allow those users
having a "rack" or parallel I/O system to take profit of that.
------------------------------------------------------------------*/

void Q2(struct params *p, struct Files* fList){
   int    *Index;
   float *dataOut;
   float *AvG; // global Average by row (now a matrx)
   int i,j;
   int nG=p->nG;
   int nE=p->nE;
   float checkVal=0;
   char nfname[500];
   FILE *fI, *ff;

   // Memory===========================================
   // Index array 
   if ((Index=(int *)calloc(nG,sizeof(int)))==NULL) 
       terror("memory for index array");
   if ((AvG =(float *)calloc(nG,sizeof(float)))==NULL)
       terror("memory for Average array");

   if ((dataOut =(float *)calloc(nG,sizeof(float)))==NULL)
       terror("memory for dataOut array");

   // Load Average vector------------
   if (p->Verbose) fprintf(stderr,"[Load Average Vector]");
   sprintf(nfname, "%s.avg",p->fOutName);
   if ((ff=fopen(nfname,"rb"))==NULL) 
        terror("opening Average Vector");
   fread(AvG, sizeof(float), nG, ff);
   fclose(ff);
   // ----------------------------------

   if ((fI=fopen("~tmp","rb"))==NULL) terror("opening tmp-index file");
   if ((ff=fopen(p->fOutName,"wb"))==NULL) terror("opening OUT file");

   // Load / compose / store================================================
   if (p->Verbose) fprintf(stderr,"[1st-Step]");

   for (i=1; i< nE; i++) { // Qnorm for each datafile: STEP 1
       fseek(fI, nG*i*sizeof(int), SEEK_SET);
       fread(Index, sizeof(int), nG, fI);
       // complete the output vector
       for (j=0;j<nG;j++) 
          dataOut[Index[j]]=AvG[j]; 

       fseek(ff, nG*i*sizeof(int), SEEK_SET);
       fwrite(dataOut, sizeof(float), nG, ff);
           
   }
   fclose(ff);
   fclose(fI);

   if (p->Verbose) fprintf(stderr,"done...\n");   

   return;

} // end Q2----------------------
 


