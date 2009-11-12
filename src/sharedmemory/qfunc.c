/******************************************************************************
* general functions for Qnormalization


* LAST REVISED: 23/02/09
******************************************************************************/
#include "qfunc.h"


// read command line arguments

struct params *command_line(int argc, char *argv[]) {
  int i;
  char c;
  struct params *parameters;
  parameters=(struct params*)malloc(sizeof(struct params));

  /* default values-------------------------------------*/
  strcpy(parameters->file_list,"qInput.txt");
  strcpy(parameters->file_out, "qOut.bin");
  parameters->transpose				= 0;
  parameters->mem_index				= 0;
  parameters->num_processors      	= DEFAULT_PROCESSOR_NUMBER;
  parameters->num_genes 			= DEFAULT_NUM_GENES;
  parameters->num_experiments		= DEFAULT_NUM_EXPERIMENTS;
  parameters->verbose  				= 0;

  /*----------------------------------------------------------*/

  for (i=1; i<argc; i++)  {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') terror("dash needed at argument (sintx error)");

    c = toupper(argv[i][1]);

    // On-Off flags----------------
    if (c == 'T') {
      parameters->transpose =1;
      continue;
    }
    if (c == 'V') {
      parameters->verbose  =1;
      continue;
    }
    if (c == 'M') {
      parameters->mem_index=1;
      continue;
    }

    // -argum=value
    if (argv[i][2]!='=') terror("equal symbol needed at argument (sintx error)");

    switch (c) {
    case 'I':
      strcpy(parameters->file_list,&argv[i][3]);
      break;
    case 'O':
      strcpy(parameters->file_out, &argv[i][3]);
      break;
    case 'E':
      parameters->num_experiments		 = atoi(&argv[i][3]);
      break;
    case 'G':
      parameters->num_genes 			 = atof(&argv[i][3]);
      break;
    case 'P':
      parameters->num_processors   	 = atoi(&argv[i][3]);
      break;
    default:
      terror("Unknown parameter (see syntax)");
    }
  }
  return parameters;
}


// panic message-----------------
void terror(char *s) {
  fprintf(stdout,"\n[ERR**] %s **\n",s);
  exit(-1);
}

// warning message---------------
void alert(char *s,char *s1) {
  fprintf(stdout,"\n[WARNING] %s : %s\n",s,s1);
}


// Load to memory a list of files
// datafile format: fileName[tab]nGenes[tab][format][newLINE]

struct files* load_list_files(struct params *parameters) {
  FILE *file_list;
  struct files*list_of_files=NULL;
  char line[MAX_SIZE_LINE],lin2[MAX_SIZE_LINE],t;
  int N=0,j,g;

  if ((file_list=fopen(parameters->file_list,"rt"))==NULL) terror("opening input file");

  if ((list_of_files=(struct files*)calloc(parameters->num_experiments,sizeof(struct files)))==NULL)
    terror("memory for list of files");

  fgets(line,MAX_SIZE_LINE,file_list);
  while (!feof(file_list)) {

    if (line[0]!='@') {

      j=sscanf(line,"%s\t%d\t%c\n",lin2,&g,&t);
      if (N==parameters->num_experiments) {
        if (parameters->verbose)
          fprintf(stderr,"[WARNING] more than %d lines in file LIST... using first %d as Experiments\n",N,N);
        parameters->num_experiments=N;
        return list_of_files;
      }

      if (lin2[strlen(lin2)-1]=='\t'||lin2[strlen(lin2)-1]==' ') lin2[strlen(lin2)-1]=0;

      if (strlen(lin2)>0) {
        list_of_files[N].fname=(char*)strdup(lin2);
        list_of_files[N].num_genes   =g;
        list_of_files[N].fType=t;
        list_of_files[N].pos=N;
        if (g > parameters->num_genes) {
          if (parameters->verbose)
            alert("more genes in file than in parameter, work with min value",list_of_files[N].fname);
          list_of_files[N].num_genes   =parameters->num_genes;
        }
        if (g < parameters->num_genes) {
          if (parameters->verbose)
            alert("more genes in parametern than in file, work with min value",list_of_files[N].fname);
          parameters->num_genes = list_of_files[N].num_genes;
        }

        N++;
      }
    }
    fgets(line,MAX_SIZE_LINE,file_list);

  }
  fclose(file_list);
  if (N!=parameters->num_experiments) {
    if (parameters->verbose) fprintf(stderr,"[WARNING] only %d files.. nExp=%d\n",N,N);
    parameters->num_experiments = N;
  }

  return list_of_files;
}



// recursive Qsort--------------------------------(not used)
void qsort_c(float *array,int l,int r,int *index) {
  int j;
  if ( l < r ) {
    // divide and conquer
    j = partition( array, l, r,index);
    //  j=(l+r)/2;
    qsort_c( array, l, j-1,index);
    qsort_c( array, j+1, r,index);
  }

}


int partition( float* a, int l, int r, int *indexes) {
  int i=l;
  int j=r+1;
  int k;
  float t,pivot;
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


// Transpose from Disk to Disk the binary matrix into a tab delimited text file

int transpose_bin_txt(struct params *parameters,char**probe_id) {
  FILE *f,*f1;
  float **mat;
  int i,j;
  char new_name[MAX_SIZE_LINE];
  int num_gen=parameters->num_genes;
  int num_experiment=parameters->num_experiments;
  int blq=1000;
  int from=0,pos;

  if (blq>num_gen) blq=num_gen;

  if ((f=fopen(parameters->file_out,"rb"))==NULL)
    terror("[Bin2Text] opening binary output file");

  sprintf(new_name,"%s.txt",parameters->file_out);
  if (parameters->verbose)fprintf(stderr,"out...%s\n",new_name);

  if ((f1=fopen(new_name,"wt"))==NULL)
    terror("[Bin2Text] opening tabulated text file out");

  if ((mat=(float**)calloc(num_experiment,sizeof(float*)))==NULL)
    terror("[Bin2Text] memory for index1");

  for (i=0; i<num_experiment;i++)
    if ((mat[i]=(float*)calloc(blq,sizeof(float)))==NULL)
      terror("[Bin2Text] memory for index2 full matrix");

  while (from<=num_gen) {
    if (parameters->verbose) fprintf(stderr,"Blq=%d of %d\n",from,num_gen);
    for (i=0; i<num_experiment;i++) { // read a BLQ of genes from each Exp
      pos=(i*num_gen+from)*sizeof(float);
      fseek(f, pos, SEEK_SET);
      fread(mat[i], blq, sizeof(float), f);
    }

    // Save a BLQ trasposed text tabulated

    for (i=0; i<blq&&from+i<num_gen;i++) {
      fprintf(f1,"%s",probe_id[from+i]);
      for (j=0; j<num_experiment; j++)
        fprintf(f1,"\t%lf",mat[j][i]);
      fprintf(f1,"\n");
    }
    from+=blq;
  }

  fclose(f);
  fclose(f1);

  return 1;
}



// Load a text-tab-2cols file ans store a bin file
int text_to_bin(struct params *parameters, char** probe_id, struct files *file_list) {

  FILE *f;
  FILE *f2;
  float*mat;
  int i;
  char new_name[MAX_SIZE_LINE];
  int num_gene=parameters->num_genes;
  int num_experiment=parameters->num_experiments;


  // probeID file
  if (parameters->verbose)fprintf(stderr,"probesID fileout...%s\n",parameters->file_out);
  if ((f=fopen(parameters->file_out,"wt"))==NULL)
    terror("[Bin2Text] opening probeID file");
  for (i=0;i<num_gene; i++) fprintf(f,"%s\n",probe_id[i]);
  fclose(f);

  if ((mat=(float *)calloc(num_gene,sizeof(float)))==NULL)
    terror("[Txt2Bin] memory for probeID array");
  if ((f2=fopen("qInBIN.txt","wt"))==NULL)
    terror("[Text2Bin] opening qInBIN file");

  if (parameters->verbose) fprintf(stderr,"kickoff..%d files.\n",num_experiment);
  for (i=0; i< num_experiment; i++) { // Qnorm for each datafile: STEP 1
    load_file(file_list, i, mat);

    if (parameters->verbose) fprintf(stderr,"%4d file: %s\n",i, file_list[i].fname);

//          sprintf(NewName,"%s.bin",fList[i].fname);
    file_list[i].fname[24]=0;
    sprintf(new_name,"../files/bin/%s.bin",&file_list[i].fname[15]);
    fprintf(f2,"%s\t6553600\tb\n",new_name);
    fflush(f2);
    if ((f=fopen(new_name,"wb"))==NULL)
      terror("[Bin2Text] opening bin file");
    fwrite(mat, sizeof(float), num_gene, f);
    fclose(f);
  }
  free(mat);
  fclose(f2);
  return 1;
}


// =======================================================================
// Dummy functions : read a datafile (case t: text: one line/one value)

void load_file(struct files*file_list, int col, float*data_in) {
  FILE *file_loaded;

  int i;
  char probe_id[255];
  float x;

  switch (file_list[col].fType) {
  case 't':
    if ((file_loaded=fopen(file_list[col].fname,"rt"))==NULL) {
      fprintf(stderr,"file:%s\n",file_list[col].fname);
      terror("opening input file (LoadFIle function)");
    }
    for (i=0;i<file_list[col].num_genes;i++) {
      fscanf(file_loaded,"%s\t%f\n",probe_id,&x);
      data_in[i]=x;
    }
    fclose(file_loaded);
    break;
  case 'b':
    if ((file_loaded=fopen(file_list[col].fname,"rb"))==NULL) {
      fprintf(stderr,"file:%s\n",file_list[col].fname);
      terror("opening input file (LoadFIle function)");
    }
    fread(data_in, file_list[col].num_genes, sizeof(float), file_loaded);
    fclose(file_loaded);
    break;
  default :
    terror("unknow (LoadFile) type of file");
  }

}

char ** load_probe_id(struct files* file_list, struct params *parameters) {
  FILE *f;
  int i;
  char probe_id[255],**probe;
  int col=0;
  float x;

  if ((probe=(char **)calloc(parameters->num_genes,sizeof(char*)))==NULL)
    terror("[Bin2Text] memory for probes");
  switch (file_list[col].fType) {
  case 't':
    if ((f=fopen(file_list[col].fname,"rt"))==NULL) {
      fprintf(stderr,"file:%s\n",file_list[col].fname);
      terror("opening input file (LoadFIle function)");
    }
    for (i=0;i<file_list[col].num_genes;i++) {
      fscanf(f,"%s\t%f\n",probe_id,&x);
      probe[i]=strdup(probe_id);
    }
    fclose(f);
    break;
  default :
    terror("unknow (loadProbeID) type of file");
  }
  return probe;

}


// DEBUG function for arrays
void debug_print(char *s, int tid,float* d, int n) {
  int j;
  fprintf(stderr, "%s----[tid=%d]--------\n",s,tid);
  for (j=0;j<n;j++) fprintf (stderr,"%f ", d[j]);
  fprintf(stderr,"\n");
}



