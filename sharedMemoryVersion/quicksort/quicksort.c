/******************************************************************************
* FILE: quicksort.c
* DESCRIPTION:
*   parallel libraries
*
* AUTHOR: J.M Mateos (13 Apr.09)
******************************************************************************/
#include "quicksort.h"

void interchange(double *x,double *y)        // swap
{
    double temp;
        
    temp=*x;
    *x=*y;
    *y=temp;
    
     
}


void interchangeIndex(int * dIndex, int index1, int index2) {
     int k = dIndex[index1];
     dIndex[index1] = dIndex[index2];
     dIndex[index2] = k;     
}


void split(double * array,int * dIndex, int first,int last,int *splitpoint)
{
    int x,i,j,s,g;
    
   
    
    // here, atleast three elements are needed
    if (array[first]<array[(first+last)/2]) {  // find median
        s=first;
        g=(first+last)/2;
    }
    else {
        g=first;
        s=(first+last)/2;
    }
    if (array[last]<=array[s]) 
        x=s;
    else if (array[last]<=array[g])
        x=last;
    else
        x=g;
    interchange(&array[x],&array[first]);      // swap the split-point element
    interchangeIndex(dIndex,x,first);   
    
                              
    
    x=array[first];
    i=first+1;                               // initialise
    j=last+1;
    while (i<j) {
        do {                                 // find j 
            j--;
        } while (array[j]>x);
        do {
            i++;                             // find i
        } while (array[i]<x);
        interchange(&array[i],&array[j]);      // swap
        interchangeIndex(dIndex,i,j);
       
       
        
    }
    
   
    interchange(&array[i],&array[j]);          // undo the extra swap
    interchangeIndex(dIndex,i,j);
    
    interchange(&array[first],&array[j]);      // bring the split-point 
    interchangeIndex(dIndex,first,j);
         
          // element to the first
    *splitpoint=j;
   
    
}





void push(int a,int b)                        // push
{
    top++;
    s[top].a=a;
    s[top].b=b;
   
}

void pop(int *a,int *b)                       // pop
{
    *a=s[top].a;
    *b=s[top].b;
    top--;
}

void insertion_sort(double * array,int *index, int first,int last)
{
    int i, c;
    double j;
    int k;
        
    for (i=first;i<=last;i++) {
        j=array[i];
        c=i;
        while ((array[c-1]>j)&&(c>first)) {
            array[c]=array[c-1];                 
            k = index[c];
            index[c] = index[c-1];
            index[c-1] = k;
            c--;
        }
        array[c]=j;
        
    }
}

void quicksort(double * array,int * dIndex, int size)
{
    int first,last,splitpoint;
    int j;
    
    s = (struct stack *) malloc((int)sizeof(struct stack)*size);
    
    
    for (j=0;j < size;j++) dIndex[j] = j;
        
    push(0,size-1); 
    while (top!=-1) {
        pop(&first,&last);
       
        for (;;) {
            
            if (last-first>SMALLSIZE) {  
               
                // find the larger sub-list
                split(array,dIndex,first,last,&splitpoint); // 
               
                // push the smaller list
                if (last-splitpoint<splitpoint-first) {
                    push(first,splitpoint-1);
                    first=splitpoint+1;
                }
                else {
                    push(splitpoint+1,last);
                    last=splitpoint-1;
                }
            }
            else {  // sort the smaller sub-lists
                    // through insertion sort
               
                insertion_sort(array,dIndex,first,last);
                
                break;
            }
        }
    }                        // iterate for larger list
    
    free(s);
}

