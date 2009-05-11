/******************************************************************************
* FILE: Qnorm.c
* DESCRIPTION:
*   Improved sequential prototype for Qnorm 
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: Jose Manuel Mateos (23 Feb.09)
* 23.Feb.09  : Using command line argums
*              Qnorm fMatrix.IN nRow nExp Normalised.fname  mode
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
* LAST REVISED: 27/02/09
******************************************************************************/

#include "Qfunc.c"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quicksort/quicksort.c"
#include "time.h"

int BLOCK_SIZE = 4;

int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;
        int myID;
        int nProcesors;
        printf("Init program \n");
        p = CommandLine(ac,av);
        
        
        MPI_Init(&ac,&av);
  
        MPI_Comm_rank(MPI_COMM_WORLD, &myID);
  
        MPI_Comm_size(MPI_COMM_WORLD,&nProcesors);

	    if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");
           
        if (myID == 0) { // Call master
             master(p,fList,nProcesors-1);    
        } else { //Call slave
           slave(p,fList,myID);  
        }
	
	return 0;
}

// Master
int master(struct params *p, struct Files* fList,int nProcesors) {
     
   int nG=p->nG;
   int nE=p->nE;
   int i;
   int index1, index2;
   int count = 0;
   int nESend = nE / nProcesors;
   int tbuf;
   char * buf;
   MPI_Status status; 
   struct Average *totalAvG; 
   struct Average *parcialAvG; 
   FILE * fI;
   FILE * fOut;
   double * dataOut;
   int numBlocks = 0;
   int limit = 0;
   int numBlocksSend =0;
   time_t start, stop;
   index1 = 0;
   index2 = 0;
   
   //time(&start);
   int index[2];
   
   float numBlocksF = 0;  
   
   if (BLOCK_SIZE > nE) {                  
      numBlocks = 1;               
   } else{
      
      numBlocksF = (float) nE / (float) BLOCK_SIZE;
   }
   
   numBlocks = (int) ceil(numBlocksF);
  // printf("nE=%i Block Real %i \n",nE,numBlocks);

   if (nProcesors  > numBlocks) {                  
       limit = numBlocks;           
   } else {          
       limit = nProcesors;   
   }
   
   for (i = 0; i < limit;i++) {
     
        index2 = index1 + BLOCK_SIZE -1 ;  
        
	if (index2 > (nE -1)) {
            index2 = nE -1;     
        } 
        
        index[0] = index1;
        index[1] = index2;
        
       // printf("Send processor %i , index1=%i index2=%i \n",i+1,index1,index2);
          
        MPI_Send(index,2,MPI_INT,i+1,1,MPI_COMM_WORLD);
        numBlocksSend++; 
        index1 = index1 + BLOCK_SIZE;      
   }
  
    
   tbuf = nG * sizeof(struct Average);
        
   buf = (char *) malloc(tbuf);   
   
   if ((totalAvG   =(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) terror("memory for Average array");  
   
   for (i=0; i <nG;i++) {
       totalAvG[i].Av=0;        // =HUGE_VAL; ???
       totalAvG[i].num=0;       
   }  
   
   //time(&stop);

   //printf("Time init slave and send initial blocks = %.0f secs \n",difftime(stop,start));
   
   //time(&start);
   while (count != numBlocks) {
         
       MPI_Recv(buf,tbuf,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
      
        
       //printf("Block slave %i \n",status.MPI_SOURCE);     
         
       parcialAvG = (struct Average *)buf;
       
       for (i = 0; i < nG; i++) {        
         totalAvG[i].Av += parcialAvG[i].Av;
         totalAvG[i].num++;           
       }      
       count++;   
       //printf("Block Recv= %i, numBlocks= %i \n",count,numBlocks);
       if (numBlocksSend < numBlocks) {
                
                index2 = index1 + BLOCK_SIZE -1 ;  
                
                if (index2 > (nE -1)) {
                    index2 = nE -1;     
                } 
                index[0] = index1;
                index[1] = index2;
      
                MPI_Send(index,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);                   
       }
         
   }
   //time(&stop);
  //printf("Time all blocks=%.0f \n",difftime(stop,start)); 
   //Send final signal all nodes
   for (i=1;i <=nProcesors;i++) {       
       index[0]= -1;
       index[1]= -1;
       MPI_Send(index,2,MPI_INT,i,1,MPI_COMM_WORLD);        
   }
   
   // Row average  ----------------------------------------------   
   for (i=0;i<nG;i++) {
       totalAvG[i].Av /=totalAvG[i].num;
      // printf("Final value %f \n",  totalAvG[i].Av);
   } 
  

   //Finally produce the ORDERED output file [SETP 2] -------------------------
   
   if ((fI=fopen("index.tmp","rb"))==NULL) terror("opening OUTPUT file");
   
   if ((fOut = fopen(p->fOutName,"wb"))==NULL) terror("opening OUTPUT file");

   if ((dataOut=(double *) calloc(nG,sizeof(double)))==NULL) terror("memory for dataout array");
    int dIndex[nG];
    int j;
    for (i = 0; i <nE;i++) {
	fseek(fI,nG*i*sizeof(int),SEEK_SET);
       	fread(dIndex,sizeof(int),nG,fI);
        //printf("Leo el fichero temp \n"); 
   	for (j=0;j<nG;j++) {
		dataOut[dIndex[j]] = totalAvG[j].Av;
	}

	fseek(fOut,(long)nG*i*sizeof(double),SEEK_SET);
	fwrite(dataOut,sizeof(double),nG,fOut);
    
    }

   fclose(fOut);
   fclose(fI);
   
   //TransposeBin2Txt(p);
  
    return 0; 
     
}


int slave(struct params *p, struct Files* fList, int myID){

   double *dataIn, *dataOut;
   int **mIndex;
   int *dIndex;
   struct Average *AvG; // global Average by row
   int i,j,k;
   FILE *fI, *fOut;
   int nG=p->nG;
   int nE=p->nE;
   int index1, index2;
   FILE *tempFile; 
   char nameFile[20];
   int fin = 1;
   time_t start, stop;
   time_t startRead, stopRead, startQuick, stopQuick; 
     MPI_Status status;  
     
     int index[2];
     
     while (fin != 0) {
     
      MPI_Recv(index,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
     // time(&start);
     // printf("Slave procesor id=%i index1=%i index2=%i \n",myID,index[0], index[1]);
     
     index1 = index[0];
     index2 = index[1];
     if ((index1 == -1) && (index2 == -1)) {
          fin = 0;       
     } else {
     
            nE = (index2 - index1) + 1;
     
            sprintf(nameFile,"index.tmp",myID);
            
            if ((tempFile=fopen(nameFile,"wb"))==NULL) terror("opening temp-index file");    
     
             // This will always be necessary to decuple the function
             if((dIndex=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

             if ((dataIn=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataIn array");
             
             if ((AvG   =(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) terror("memory for Average array");
   
             for (j=0; j< nG;j++) { // init Accumulation array 
                  AvG[j].Av=0;        // =HUGE_VAL; ???
                  AvG[j].num=0;
             }

       
               // QNORM ===============================================================
               for (i=0; i< nE; i++) { // Qnorm for each datafile: STEP 1
                    // time(&startRead);
		     LoadFile(fList, i+index1, dataIn,nG); //Load a file
                   // time(&stopRead);
            #ifdef DEBUG
                   printf("Show file ");
                   DebugPrint("Load", dataIn, fList[i+index1].nG); 
            #endif
                   // time(&startQuick);
                    Qnorm1(dataIn, dIndex, nG); // dataIn returns ordered and Index contains the origial position
                   // time(&stopQuick);
                   
		    #ifdef DEBUG
                    DebugPrint("Sorted", dataIn, nG);         
            #endif
                    
                    AccumulateRow(AvG, dataIn , nG); // Calcula la media
            
                    
                   // for (j=0;j<nG;j++) {
                   //     mIndex[j][i]= dIndex[j]; // Aqui tengo que cambiar la variable i
                   // }
            	fseek(tempFile,(i+index1)*nG*sizeof(int),SEEK_SET);
            	fwrite(dIndex,sizeof(int),nG,tempFile);
            
            #ifdef DEBUG
                   // fprintf(stderr,"Index (col=%d)\n",i);
                   // for (j=0;j<nG;j++) fprintf (stderr,"%d ", dIndex[j]); fprintf(stderr,"\n");
            #endif
            
                  //  printf("Times myID=%i file=%i time Read=%.0f time Quick=%.0f \n",myID,index1+i,difftime(stopRead,startRead),difftime(stopQuick,startQuick));
            
               } // End bucle
   

                fclose(tempFile);
   
               // Row average  ----------------------------------------------   
               for (i=0;i<nG;i++) {
                   AvG[i].Av /=AvG[i].num;
               }             
               
               char *buf;
               
               if ((buf = (char *) malloc(sizeof(struct Average)* nG)) == NULL)  terror("ERROR MEMORY: Sending diagonal");   
               buf = (char *) AvG;
               time(&stop);
               printf("Send data from slave id=%i tiempo=%.0f   \n",myID,difftime(stop,start));
               MPI_Send(buf,sizeof(struct Average)*nG,MPI_BYTE,0,1,MPI_COMM_WORLD);
               
     }  
   }

   

   return 0;

}


// input returns ordered and Index contains the origial position

int Qnorm1(double *input, int *dIndex, int nG){
	int i,j,k,n;
       time_t start, stop; 
       // time(&start);
	for (j=0; j<nG;j++) dIndex[j]=j; // init the indexes array

/*
	for (j=0; j<nG;j++) // UNIFY NAN CONSTANT 
	   if ((!__finite(input[j]))||(__isnan(input[j]))) input[j]=HUGE_VAL;
*/

//	QsortC(input,0,nG-1,dIndex); // Quicksort 
	
        
        
	quicksort(input,dIndex,nG);
       // time(&stop);
       // printf("Time quick en function = %.0f \n",difftime(stop,start));
        return 1;
}

void AccumulateRow(struct Average *AvG, double *input , int nG){
        int i;
 	
        for (i=0;i<nG;i++) {

/*
           if ((__finite(input[i])&&(!__isnan(input[i])))){
		if ((!__finite(AvG[i].Av))||(__isnan(AvG[i].Av))){
			AvG[i].Av=0;
		}
*/
		AvG[i].Av+=input[i];
		AvG[i].num++;
/*           }
*/	}


	return;

}


