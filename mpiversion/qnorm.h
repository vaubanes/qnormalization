/*
 * qnorm.h
 *
 *  Created on: 11/08/2009
 *      Author: bitlab
 */

#ifndef QNORM_H_
#define QNORM_H_

#pragma once

#include "qfunc.h"
#include "qsort.h"

int main(int ac, char **av);
int master(struct params *p, struct Files* fList,int nProcesors);
int slave(struct params *p, struct Files* fList, int myID);



#endif /* QNORM_H_ */
