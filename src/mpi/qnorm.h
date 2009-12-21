/*
 *FILE:qnorm.h
 *
 *
 *  Author: Jose Manuel Mateos
 */

#ifndef QNORM_H_
#define QNORM_H_


#include "qfunc.h"


int main(int ac, char **av);
int master(Params *p,int nProcesors);
int slave(Params *p, InfoFile* fList, int myID);
int store_final_result(Average *total_avg, int num_experiments,int num_genes, char * name_out_file);



#endif /* QNORM_H_ */
