CC = gcc
CFLAGS = -Wall -c -fopenmp
PFLAGS = -lgomp
OBJS = qnorm_a.o qfunc.o quicksort.o

qnorm_a.o : qnorm_a.c
	${CC} ${CFLAGS} qnorm_a.c

qfunc.o: qfunc.c
	${CC} ${CFLAGS} qfunc.c

quicksort.o: quicksort.c
	${CC} ${CFLAGS} quicksort.c

shared_qnorm: ${OBJS}
	${CC} ${PFLAGS} ${OBJS} -o shared_qnorm

all: shared_qnorm

defaul_target: all

clean:
	rm -f ${OBJS}

test: shared_qnorm
	./shared_qnorm -i=qInBIN.txt -e=4 -g=80 -o=qOutBin.bin -P=1

  




