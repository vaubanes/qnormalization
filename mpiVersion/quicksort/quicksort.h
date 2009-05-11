/******************************************************************************
* FILE: quicksort.h
* jmateos@uma.es
* LAST REVISED: 13/04/09
******************************************************************************/

#include <stdio.h>
#include "time.h"


#define SMALLSIZE       100            // not less than 3


struct stack {                              // stack element.
        int a,b;
} * s;

int top=-1;  
