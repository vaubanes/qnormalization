/*
 * qnorm_a.h
 *
 *  Created on: 16/10/2009
 *      Author: bitlab
 */

#ifndef QNORM_A_H_
#define QNORM_A_H_

#include "qfunc.h"
#include "quicksort.h"

int main(int ac, char **av);
void qnorm_main(struct params *p, struct Files* fList);
int qnorm_sort(float*input, int *dIndex, int nG);
void accumulate_row(struct Average *AvG, float*input , int nG);





#endif /* QNORM_A_H_ */
