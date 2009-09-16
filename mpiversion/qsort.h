/*
 * qsort.h
 *
 *
 *      Author: Jose Manuel Mateos
 */

#ifndef QSORT_H_
#define QSORT_H_



#define SMALLSIZE       10            // not less than 3

typedef struct {                              // stack element.
  int a,b;
}Stack;


void interchange(double *x, double *y);
void interchange_index(int * dIndex, int index1,int index2);
void split(double * array,int * dIndex, int first,int last,int *splitpoint);
void push(int a,int b);
void pop(int *a,int *b);
void insertion_sort(double * array,int *index, int first,int last);
void quicksort(double * array,int * dIndex, int size);

#endif /* QSORT_H_ */
