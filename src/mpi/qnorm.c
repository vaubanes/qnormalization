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
*  -n      Block Size 				   2               positive integer
* ---------------------------------------------------------------------------
*
*	@returns >0 if everything was fine <0 if there was an error
*
* LAST REVISED: 14/08/09
******************************************************************************/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "time.h"
#include "qnorm.h"


int main(int ac, char **av) {

  InfoFile *flist=NULL;
  Params *parameters=NULL;
  int myid;
  int num_processors;

  MPI_Init(&ac,&av);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  MPI_Comm_size(MPI_COMM_WORLD,&num_processors);

  printf("Processor ID=%i working Processors=%i \n ",myid,num_processors);

  parameters = commandline(ac,av);

  if ((flist=load_input_files(parameters))==NULL)
    terror("Loading list of files");

  if (myid == 0) { // Call master
    master(parameters,num_processors-1);
  } else { //Call slave
    slave(parameters,flist,myid);
  }

  printf("Processor ID=%i finish \n",myid);
  MPI_Finalize();
  return 0;
}



// Execute in master processor
int master(Params *p,int num_processors) {

  const int num_genes=p->num_genes;
  const int num_experiments=p->num_experiments;
  const int block_size = p->block_size;
  int i;
  int index1, index2;
  int count = 0;
  int text_buf;
  char * buf;
  MPI_Status status;
  Average *total_avg;
  Average *parcial_avg;
  FILE * file_input;
  FILE * file_out;
  double * dataout;
  int numblocks = 0;
  int limit = 0;
  int numblocks_send =0;
  char namefile[20];
  int index[2];

  index1 = 0;
  index2 = 0;
  // Get number of blocks
  numblocks =  calculate_blocks(num_experiments,block_size);

  limit =calculate_initial_blocks(numblocks, num_processors);

  index[0] = -1;
  index[1] = -1;

  // Send initials blocks
  for (i = 0; i < limit;i++) {

    calculate_index_blocks(index,num_experiments,block_size);

    MPI_Send(index,2,MPI_INT,i+1,1,MPI_COMM_WORLD);

    numblocks_send++;
  }



  text_buf = num_genes * sizeof(Average);

  if ((buf = (char *) malloc(text_buf))==NULL) terror("memory buffer array");

  if ((total_avg   =(Average *)calloc(num_genes,sizeof(Average)))==NULL) terror("memory for Average array");

  // Inicialize average array
  for (i=0; i <num_genes;i++) {
    total_avg[i].value=0;
    total_avg[i].elements=0;
  }


  while (count != numblocks) {

    MPI_Recv(buf,text_buf,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

    parcial_avg = (Average *)buf;

    // Acumulate partial average
    for (i = 0; i < num_genes; i++) {
      total_avg[i].value += parcial_avg[i].value;
      total_avg[i].elements+= parcial_avg[i].elements;
    }

    count++;

    if (numblocks_send < numblocks) {

      calculate_index_blocks(index,num_experiments,block_size);
      numblocks_send++;
      MPI_Send(index,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
    }

  }

  //Send final signal all nodes
  for (i=1;i <=num_processors;i++) {
    index[0]= -1;
    index[1]= -1;
    MPI_Send(index,2,MPI_INT,i,1,MPI_COMM_WORLD);
  }

  // Calculate final  average  ----------------------------------------------
  for (i=0;i<num_genes;i++) {
    total_avg[i].value /=total_avg[i].elements;
  }



  if ((file_out = fopen(p->file_out,"wb"))==NULL) terror("opening OUTPUT file");

  if ((dataout=(double *) calloc(num_genes,sizeof(double)))==NULL) terror("memory for dataout array");

  int data_index[num_genes];
  int j;

  char command[1024];

  for (i = 0; i <num_experiments;i++) {

    sprintf(namefile,"index%i.tmp",i);
    if ((file_input=fopen(namefile,"rb"))==NULL) terror("opening Index file");

    fseek(file_input,0,SEEK_SET);
    fread(data_index,sizeof(int),num_genes,file_input);
    fclose(file_input);

    sprintf(command,"rm index%i.tmp",i);
    system(command);


    for (j=0;j<num_genes;j++) {

      dataout[data_index[j]] = total_avg[j].value;
    }

    fseek(file_out,(long)num_genes*i*sizeof(double),SEEK_SET);
    fwrite(dataout,sizeof(double),num_genes,file_out);

  }

  fclose(file_out);

  if (p->transpose) {

    transpose_matrix(p);
  }

  return 0;

}


int slave(Params *p, InfoFile* flist, int myid) {

  double * data_input;
  int *dIndex;
  Average *average; // global Average by row
  int i,j;

  int num_genes=p->num_genes;
  int num_experiments=p->num_experiments;
  int index1, index2;
  FILE *index_out;
  char name_file[20];
  int end = 1;
  MPI_Status status;
  int index[2];

  while (end != 0) {
    //Recieve index of experiments from master
    MPI_Recv(index,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status); //

#if DEBUG
    // printf("Slave procesor id=%i index1=%i index2=%i \n",myID,index[0], index[1]);
#endif

    index1 = index[0];
    index2 = index[1];

    if ((index1 == -1) && (index2 == -1)) { //if value index is -1 mind master send end signal to slave
      end = 0;
    } else {

      num_experiments = (index2 - index1) + 1; // Calculate number of experiments to work

      if ((dIndex=(int *)calloc(num_genes,sizeof(int)))==NULL) terror("memory for index2");

      if ((data_input=(double *)calloc(num_genes,sizeof(double)))==NULL) terror("memory for dataIn array");

      if ((average=(Average *)calloc(num_genes,sizeof(Average)))==NULL) terror("memory for Average array");

      for (j=0; j< num_genes;j++) { // init Accumulation array
        average[j].value=0;
        average[j].elements=0;
      }


      for (i=0; i< num_experiments; i++) { // Apply Qnorm for each experiments

        load_parcial_result(flist, i+index1, data_input,num_genes); //Load value of genes of experiment i

#ifdef DEBUG
        printf("Show file ");
        debug_print("Load", data_input, flist[i+index1].num_genes);
#endif

        qnorm(data_input, dIndex, num_genes); // data_input returns ordered and Index contains the origial position

#ifdef DEBUG
        debug_print("Sorted", data_input, num_genes);
#endif

        accumulate_row(average, data_input , num_genes); // Calulate partial average

        // Save result in a temporal file
        sprintf(name_file,"index%i.tmp",i+index1);
        if ((index_out = fopen(name_file,"wb"))==NULL) terror("ERROR: Open file");

        //Write array with index array with original positions;
        fwrite(dIndex,sizeof(int),num_genes,index_out);

        fclose(index_out);
      } // End  bucle for

      char *buf;

      if ((buf = (char *) malloc(sizeof(Average)* num_genes)) == NULL)  terror("ERROR MEMORY: Sending diagonal");
      buf = (char *) average;

      MPI_Send(buf,sizeof(Average)*num_genes,MPI_BYTE,0,1,MPI_COMM_WORLD);

    }

  }

  return 0;

}




