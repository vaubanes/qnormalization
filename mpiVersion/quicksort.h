/******************************************************************************
* FILE: quicksort.h
* jmateos@uma.es
* LAST REVISED: 13/04/09
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>



#define SMALLSIZE       10            // not less than 3


struct stack {                              // stack element.
  int a,b;
} * s;




void interchange(double *x, double *y);
void interchangeIndex(int * dIndex, int index1,int index2);

