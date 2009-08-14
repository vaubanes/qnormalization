/*
 * qnorm.h
 *
 *
 *      Author: Jose Manuel Mateos
 */

#ifndef QNORM_H_
#define QNORM_H_


#include "qfunc.h"


int main(int ac, char **av);
int master(struct Params *p,int nProcesors);
int slave(struct Params *p, struct Files* fList, int myID);



#endif /* QNORM_H_ */
