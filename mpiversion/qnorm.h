/*
 * qnorm.h
 *
 *
 *      Author: Jose Manuel Mateos
 */

#ifndef QNORM_H_
#define QNORM_H_


#include "qfunc.h"
#include "qsort.h"

int main(int ac, char **av);
int master(struct params *p,int nProcesors);
int slave(struct params *p, struct Files* fList, int myID);



#endif /* QNORM_H_ */
