/**
* GPU QNORM implementation
* -------------------------
**/

#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gpuqsort.h"


#include "timer.cpp"
#include "Qfunc.h"   // funciones OTS


typedef unsigned int element;  // --->  definimos el tipo del template


/**
* The program entrance point
*
* Give no argument to show information
*/
int main(int argc, char* argv[])
{

/**
*/
	// Allocate memory for the sequences to be sorted
        andresTimer st;
 	double t1,t2;
        float dato;
        unsigned int maxsize=7000000;
 	st.start();  // tiempo de carga
// OTS qnorm:
        struct Files *fList=NULL;
        struct params *p=NULL;
        p = CommandLine(argc,argv);
	if ((fList=LoadListOfFiles(p))==NULL) terror("Loading list of files");
        //QNormMain(p,fList);	
// end OTS
        
	unsigned int i=0;
        element* data = new element[maxsize];
        unsigned int* pos = new unsigned int[maxsize];
        float * cum = new float[maxsize];

printf("\n\nQnorm GPU version (%s)\n",argv[0]);


	double timerValue;

				// Store copy of sequence
				//memcpy(data2,data,testsize*sizeof(element));

				// Sort it
				if(gpuqsort(fList,p,p->blocks,p->threads,p->shared,0)!=0)
				{
					printf("Error! (%s)\n",getGPUSortErrorStr());
					exit(1);
				}

           
}
