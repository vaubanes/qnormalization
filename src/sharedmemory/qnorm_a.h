/*
 * qnorm_a.h
 *
 *  Created on: 23/02/2009
 *      Author: ots
 */

#ifndef QNORM_A_H_
#define QNORM_A_H_

#include "qfunc.h"
#include "quicksort.h"

int main(int ac, char **av);
void qnorm_main(struct params *parameters, struct files* file_list);
int qnorm_sort(float*input, int *dIndex, int nG);
void accumulate_row(struct average *average, float*input , int num_genes);
int omp_get_thread_num();
int omp_get_num_threads();





#endif /* QNORM_A_H_ */
