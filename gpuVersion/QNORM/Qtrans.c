/******************************************************************************
* FILE: pQtras.c
* DESCRIPTION: Transpose from Disk to Disk the binary matrix into a 
               tab delimited text file
               Include probeID from first file

* 23.Feb.09  : Using command line argums
*              pQtrans [-o=Value]  see below Command line params
*              (same as pQnormX.c)

*              fList  : contains a list of nExp filenames (gene-expression data files) 
*                          line format: fileName[TAB]nGenes[TAB]FileType[NEWLINE]
*              nRows  : number of genes in each file 
*              Normalised.fname where the normalised values will be stored 
*                       (as a text-tabulated matrix)
*
*  Option  Description              Default value     alternative Values
*  ------  -----------------------  -------------     ------------------
*  -i      File name (list of files)   qInput.txt      valid existing pathname
*  -o      Output binary matrix        qOut.bin        binary by columns file
*  -g      NUmber of genes             15              positive integer
*  -e      Number of experiments       2               positive integer
*  -V      Verbose mode                Not             -V  
* ---------------------------------------------------------------------------
*
* LAST REVISED: 30/03/09
******************************************************************************/
#include "Qfunc.c"

int main(int ac, char **av){

        struct Files *fList=NULL;
        struct params *p=NULL;
        char **probeID;
        fprintf(stderr,"\n[transposing Disk->Disk]\n"); fflush(stderr);
        p = CommandLine(ac,av);

	if ((fList=LoadListOfFiles(p))==NULL) 
           terror("Loading list of files");
        probeID=LoadprobeID(fList, p);

        TransposeBin2Txt(p,probeID);
	
	return 1;
}



