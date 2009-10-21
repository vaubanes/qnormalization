/******************************************************************************
* FILE: quicksort.h
* jmateos@uma.es
* LAST REVISED: 13/04/09
******************************************************************************/

#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <stdio.h>
#include <stdlib.h>



#define SMALLSIZE       10            // not less than 3


struct stack {                              // stack element.
        int a,b;
} * s;



void interchange(float *x,float *y);
void interchangeIndex(int * dIndex, int index1, int index2);
void split(float * array,int * dIndex, int first,int last,int *splitpoint);
void push(int a,int b);
void pop(int *a,int *b);
void insertion_sort(float* array,int *index, int first,int last);
void quicksort(float * array,int * dIndex, int size);

#endif
