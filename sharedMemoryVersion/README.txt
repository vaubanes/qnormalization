*====================================================================================
* FILE: Qtrans.c (v.6 : Qtrans6.c)

*       gcc Qtrans6.c -o pqtrans


* DESCRIPTION: Transpose from Disk to Disk the binary matrix into a 
               tab delimited text file
               Include probeID from first file

* 13.Abr.09  : Using command line argums (same as pQnorm)
*
*              qtrans [-o=Value]  see below Command line params
*              (same as pQnorm6.c)

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



*====================================================================================
* FILE: pQnorm.c (version 6) (pQnorm6.c -> pqnorm)

*       gcc Qnorm6.c -o qnorm

* DESCRIPTION:
*   parallel version of Quantile normalization procedure
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
******************************************************************************/


Qfunc.c : general functions for Qnormalization
Qfunc.h : headers and structures

quicksort/quicksort.c : iterative quicksort (to solve stack problems)

ejeQn: script to launch sequential Qnorm program

ejeQtrans: script to launch Transpossing matrix (disk to disk)


qIn470.txt : list of files to be proceesed (470 files)