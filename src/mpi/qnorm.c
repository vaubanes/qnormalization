/******************************************************************************
* FILE: Qnorm.c
* DESCRIPTION:
*
*   MPI version  prototype for Qnorm
*   Qnormalisation Method: function that implements the ben Bolstad Method
*   quantile Normalization of High density Oliglonucleotide Array Data
*
* AUTHOR: Jose Manuel Mateos Duran (23 Feb.09)
*
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
* LAST REVISED: 16/12/09
******************************************************************************/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "time.h"
#include "qnorm.h"

// Point init program, all processors start here
int main(int ac, char **av) {

  InfoFile *file_list_of_files=NULL;
  Params *init_parameters=NULL;
  int myid;
  int num_processors;

  // Routines to initialize mpi
  MPI_Init(&ac,&av);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD,&num_processors);

  printf("Processor ID=%i working Processors=%i \n ",myid,num_processors);

  // Read initial parameters
  init_parameters = commandline(ac,av);

  // Read file that contains a list of files of experiments
  if ((file_list_of_files=load_input_files(init_parameters))==NULL)
    terror("Loading list of files");

  // Distinction between master processor and slave processors
  if (myid == 0) { // Call master
    master(init_parameters,num_processors-1);
  } else { //Call slave
    slave(init_parameters,file_list_of_files,myid);
  }

  printf("Processor ID=%i finish \n",myid);

  MPI_Finalize();

  return 0;
}


// This function is executed for one processor
// This function coordinate the work between slave processors
int master(Params *p,int num_processors) {

  const int num_genes=p->num_genes;
  const int num_experiments=p->num_experiments;
  const int block_size = p->block_size;
  int i;
  int count = 0;
  int size_buffer;
  char * buffer_recv;
  MPI_Status status;
  Average *total_avg;
  Average *parcial_avg;
  int num_blocks = 0;
  int number_of_initial_blocks = 0;
  int numblocks_send =0;
  int index[2];

  // The master sends blocks of work to slave processors
  // these blocks are composed of a variable set of experiments

  // First we have to calculate the total number of blocks (depend on block_size)
  num_blocks =  calculate_blocks(num_experiments,block_size);

  // The master makes an initial distribution of blocks between the slaves
  // For to do this, the master needs to know the number of initial blocks to distribute
  number_of_initial_blocks =calculate_initial_blocks(num_blocks, num_processors);

  // To indicate a block, the master use a array of two elements where:
  index[0] = -1; // Indicate the first experiment of the block
  index[1] = -1; // Indicate the last experiment of the block

  // Send initials blocks to slave processors
  for (i = 0; i < number_of_initial_blocks;i++) {
	// Calculate the index of the block
    calculate_index_blocks(index,num_experiments,block_size);

    MPI_Send(index,2,MPI_INT,i+1,1,MPI_COMM_WORLD);

    numblocks_send++;
  }

  //Calculate the size of buffer
  size_buffer = num_genes * sizeof(Average);

  // The master used this buffer to recieve the result of slave
  if ((buffer_recv = (char *) malloc(size_buffer))==NULL) terror("memory buffer array");

  // This vector will contain the final averages for each gene contains number of experiments and summation
  if ((total_avg   =(Average *)calloc(num_genes,sizeof(Average)))==NULL) terror("memory for Average array");

  // Inicialize average array
  for (i=0; i <num_genes;i++) {
    total_avg[i].summation=0;
    total_avg[i].num_experiments=0;
  }

  // This while will be alive while master will not send all blocks
  while (count != num_blocks) {

	// Recieve the results
    MPI_Recv(buffer_recv,size_buffer,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

    // The buffer contains the partial average of a set of experiments
    parcial_avg = (Average *)buffer_recv;

    // Add the partial average to global average
    for (i = 0; i < num_genes; i++) {
      total_avg[i].summation += parcial_avg[i].summation;
      total_avg[i].num_experiments+= parcial_avg[i].num_experiments;
    }

    count++;

    // if master still has not sent all blocks, send the next block to last slave processor
    if (numblocks_send < num_blocks) {

      calculate_index_blocks(index,num_experiments,block_size);
      numblocks_send++;
      MPI_Send(index,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
    }

  } // end while


  //Send a end signal all slave processors
  for (i=1;i <=num_processors;i++) {
    index[0]= -1;
    index[1]= -1;
    MPI_Send(index,2,MPI_INT,i,1,MPI_COMM_WORLD);
  }

  // Calculate final  average for all genes
  for (i=0;i<num_genes;i++) {
    total_avg[i].summation /=total_avg[i].num_experiments;
  }

  // Store the final result in disk
  store_final_result(total_avg,num_experiments,num_genes,p->file_out);

  // If In the arguments is active the transpose parameter
  if (p->transpose) {
	 // Transpose the result matrix to origin positions
    transpose_matrix(p);
  }

  return 0;

}

// This function is executed for all slave processors
int slave(Params *p, InfoFile* flist, int myid) {

  double * experiment;
  int *initial_index_experiment;
  Average *block_average; // global Average by row
  int i,j;

  int num_genes=p->num_genes;
  int num_experiments=p->num_experiments;
  int index1, index2;
  FILE *index_out;
  char name_file[20];
  int end = 1;
  MPI_Status status;
  int block_indexs[2];

  // While the master doesn't send end signal
  while (end != 0) {
    //Receive from master a block, by a matrix with start and end indexs
    MPI_Recv(block_indexs,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status); //

#if DEBUG
    // printf("Slave processor id=%i index1=%i index2=%i \n",myID,index[0], index[1]);
#endif

    index1 = block_indexs[0];
    index2 = block_indexs[1];

    //if value index is -1 mind master send end signal to slave and end the while
    if ((index1 == -1) && (index2 == -1)) {
      end = 0;
    } else {

      // Calculate the number of experiment of the block
      num_experiments = (index2 - index1) + 1;

      // Reserve memory to store original index of a experiment
      if ((initial_index_experiment=(int *)calloc(num_genes,sizeof(int)))==NULL) terror("memory for index2");

      // Reserve memory to store a experiment
      if ((experiment=(double *)calloc(num_genes,sizeof(double)))==NULL) terror("memory for dataIn array");

      // Reserve memory to store average of a experiment
      if ((block_average=(Average *)calloc(num_genes,sizeof(Average)))==NULL) terror("memory for Average array");

      for (j=0; j< num_genes;j++) { // init Accumulation array
        block_average[j].summation=0;
        block_average[j].num_experiments=0;
      }

      // For each experiment into a block
      for (i=0; i< num_experiments; i++) {
    	// Load a experiment form file
        load_experiment_from_file(flist, i+index1, experiment,num_genes); //Load value of genes of experiment i

#ifdef DEBUG
        printf("Show file ");
        debug_print("Load", data_input, flist[i+index1].num_genes);
#endif
        // Sort the experiment vector and store in initial index vector the original positions
        qnorm(experiment, initial_index_experiment, num_genes);

#ifdef DEBUG
        debug_print("Sorted", data_input, num_genes);
#endif
        // Add the vector sorted to a experiment average vector
        accumulate_row(block_average, experiment , num_genes);

        // Save the result in a temporal file
        sprintf(name_file,"index%i.tmp",i+index1);
        if ((index_out = fopen(name_file,"wb"))==NULL) terror("ERROR: Open file");
        //Write array with index array with original positions;
        fwrite(initial_index_experiment,sizeof(int),num_genes,index_out);

        fclose(index_out);
      } // End  bucle for

      // When the slave has finished all experiments of the a block sends the average to master
      char *buf;

      if ((buf = (char *) malloc(sizeof(Average)* num_genes)) == NULL)  terror("ERROR MEMORY: Sending diagonal");
      buf = (char *) block_average;

      MPI_Send(buf,sizeof(Average)*num_genes,MPI_BYTE,0,1,MPI_COMM_WORLD);

    }

  }

  return 0;

}

// This function stores the results in a file.
int store_final_result(Average *total_avg, int num_experiments,int num_genes, char * name_out_file) {

	 FILE * file_input;
	 FILE * file_out;
	 int initial_positions[num_genes];
	 double * experiment;
	 int i, j;
	 char command[1024];
	 char namefile[20];

	  // Open the file
	  if ((file_out = fopen(name_out_file,"wb"))==NULL) terror("opening OUTPUT file");

	  // Reserv memory to store a experiment
	  if ((experiment=(double *) calloc(num_genes,sizeof(double)))==NULL) terror("memory for dataout array");


	  // For each experiment
	  for (i = 0; i <num_experiments;i++) {

		// Read the partial result stored by a slave
	    sprintf(namefile,"index%i.tmp",i);
	    if ((file_input=fopen(namefile,"rb"))==NULL) terror("opening Index file");

	    fseek(file_input,0,SEEK_SET);
	    fread(initial_positions,sizeof(int),num_genes,file_input);
	    fclose(file_input);

	    // Remove the temporal files created by slave processor
	    sprintf(command,"rm index%i.tmp",i);
	    system(command);

	    // For each genes of a experiment
	    for (j=0;j<num_genes;j++) {
	      // Save in the original position the average value by row
	      experiment[initial_positions[j]] = total_avg[j].summation;
	    }

	    // Store the experiment
	    fseek(file_out,(long)num_genes*i*sizeof(double),SEEK_SET);
	    fwrite(experiment,sizeof(double),num_genes,file_out);

	  } // end for

	  fclose(file_out);

	  return 0;
}




