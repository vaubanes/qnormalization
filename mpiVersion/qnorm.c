/******************************************************************************
* FILE: Qnorm.c
* DESCRIPTION:
*
*   MPI version  prototype for Qnorm
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: Jose Manuel Mateos Duran (23 Feb.09)
* 23.Feb.09  : Using command line argums
*              Qnorm fMatrix.IN nRow nExp Normalised.fname  mode
*
*              fList  : contains a list of nExp filenames (gene-expression data files)
*                          line format: fileName[TAB]nGenes[TAB]FileType[NEWLINE]
*              nRows  : number of genes in each file
*              Normalised.fname where the normalised values will be stored
*                       (as a text-tabulated matrix)
*
*
*  Command line Parameters
*         sintaxis:    Qnorm [-Option=value]...
*
*  Option  Description              Default value     alternative Values
*  ------  -----------------------  -------------     ------------------
*  -i      File name (list of files)   qInput.txt      valid existing pathname
*  -o      Output binary matrix        qOut.bin        binary by columns file
*  -e      Number of experiments       2               positive integer
*  -g      NUmber of genes             15              positive integer
*  -t      Traspose the fileOut        Not             -T (yes)
*
* ---------------------------------------------------------------------------
*
*	@returns >0 if everything was fine <0 if there was an error
*
* LAST REVISED: 27/02/09
******************************************************************************/

#include "Qfunc.c"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quicksort.c"
#include "time.h"



int main(int ac, char **av) {

  struct Files *fList=NULL;
  struct params *p=NULL;
  int myID;
  int nProcesors;

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
  index1 = 0;
  index2 = 0;
  char nameFile[20];
  int index[2];

  // Get number of blocks
  numBlocks =  calculateBlocks(nE);

  limit = calculateInitialBlocks(numBlocks, nProcesors);

  index[0] = -1;
  index[1] = -1;

  // Send initials blocks
  for (i = 0; i < limit;i++) {

    calculateIndexBlocks(index,nE);

    MPI_Send(index,2,MPI_INT,i+1,1,MPI_COMM_WORLD);

    numBlocksSend++;
  }


  tbuf = nG * sizeof(struct Average);

  if ((buf = (char *) malloc(tbuf))==NULL) terror("memory buffer array");

  if ((totalAvG   =(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) terror("memory for Average array");

  // Inicialize average array
  for (i=0; i <nG;i++) {
    totalAvG[i].Av=0;
    totalAvG[i].num=0;
  }


  while (count != numBlocks) {

    MPI_Recv(buf,tbuf,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

    parcialAvG = (struct Average *)buf;

    // Acumulate partial average
    for (i = 0; i < nG; i++) {
      totalAvG[i].Av += parcialAvG[i].Av;
      totalAvG[i].num+= parcialAvG[i].num;
    }

    count++;

    if (numBlocksSend < numBlocks) {

      calculateIndexBlocks(index,nE);

      // printf("*** Send blocks index1=%i index2=%i ID=%i \n",index[0],index[1],status.MPI_SOURCE);
      numBlocksSend++;
      MPI_Send(index,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
    }

  }


  //Send final signal all nodes
  for (i=1;i <=nProcesors;i++) {
    index[0]= -1;
    index[1]= -1;
    MPI_Send(index,2,MPI_INT,i,1,MPI_COMM_WORLD);
  }

  // Calculate final  average  ----------------------------------------------
  for (i=0;i<nG;i++) {

    totalAvG[i].Av /=totalAvG[i].num;
  }



  if ((fOut = fopen(p->fOutName,"wb"))==NULL) terror("opening OUTPUT file");

  if ((dataOut=(double *) calloc(nG,sizeof(double)))==NULL) terror("memory for dataout array");

  int dIndex[nG];
  int j,h, value;
  char line[1024];
  char command[1024];

  for (i = 0; i <nE;i++) {

    sprintf(nameFile,"index%i.tmp",i);
    if ((fI=fopen(nameFile,"rb"))==NULL) terror("opening Index file");

    fseek(fI,0,SEEK_SET);
    fread(dIndex,sizeof(int),nG,fI);
    fclose(fI);
    sprintf(command,"rm index%i.tmp",i);
    system(command);

    for (j=0;j<nG;j++) {
      dataOut[dIndex[j]] = totalAvG[j].Av;
    }

    fseek(fOut,(long)nG*i*sizeof(double),SEEK_SET);
    fwrite(dataOut,sizeof(double),nG,fOut);

  }

  fclose(fOut);


  if (p->Traspose) {

    TransposeBin2Txt(p);
  }


  return 0;

}


int slave(struct params *p, struct Files* fList, int myID) {

  double *dataIn, *dataOut;
  int **mIndex;
  int *dIndex;
  struct Average *AvG; // global Average by row
  int i,j,k;
  FILE *fI, *fOut;
  int nG=p->nG;
  int nE=p->nE;
  int index1, index2;
  FILE *tempFile,*indexOut;
  char nameFile[20];
  int fin = 1;
  MPI_Status status;
  int index[2];

  while (fin != 0) {

    MPI_Recv(index,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

#if DEBUG
    // printf("Slave procesor id=%i index1=%i index2=%i \n",myID,index[0], index[1]);
#endif

    index1 = index[0];
    index2 = index[1];

    if ((index1 == -1) && (index2 == -1)) {
      fin = 0;
    } else {

      nE = (index2 - index1) + 1;

      if ((dIndex=(int *)calloc(nG,sizeof(int)))==NULL) terror("memory for index2");

      if ((dataIn=(double *)calloc(nG,sizeof(double)))==NULL) terror("memory for dataIn array");

      if ((AvG   =(struct Average *)calloc(nG,sizeof(struct Average)))==NULL) terror("memory for Average array");

      for (j=0; j< nG;j++) { // init Accumulation array
        AvG[j].Av=0;
        AvG[j].num=0;
      }


      for (i=0; i< nE; i++) { // Qnorm for each datafile: STEP 1

        LoadFile(fList, i+index1, dataIn,nG); //Load a file

#ifdef DEBUG
        printf("Show file ");
        DebugPrint("Load", dataIn, fList[i+index1].nG);
#endif

        Qnorm1(dataIn, dIndex, nG); // dataIn returns ordered and Index contains the origial position

#ifdef DEBUG
        DebugPrint("Sorted", dataIn, nG);
#endif

        AccumulateRow(AvG, dataIn , nG); // Calulate partial average

        // Save file index
        sprintf(nameFile,"index%i.tmp",i+index1);
        if ((indexOut = fopen(nameFile,"wb"))==NULL) terror("ERROR: Open file");
        /**
        int h = 0;
        for(h=0;h <nG;h++) {

                     fprintf(indexOut,"%i \n",dIndex[h]);
        }**/

        fseek(indexOut,0,SEEK_SET);
        fwrite(dIndex,sizeof(int),nG,indexOut);

      } // End  bucle

      fclose(indexOut);


      char *buf;

      if ((buf = (char *) malloc(sizeof(struct Average)* nG)) == NULL)  terror("ERROR MEMORY: Sending diagonal");
      buf = (char *) AvG;

      MPI_Send(buf,sizeof(struct Average)*nG,MPI_BYTE,0,1,MPI_COMM_WORLD);

    }
  }


  return 0;

}




