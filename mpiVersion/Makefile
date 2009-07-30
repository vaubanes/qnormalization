# CC = gcc
CC = mpicc 
CFLAGS = -Wall
LFLAGS = -lm

OBJS = qnorm.o qfunc.o quicksort.o

qnorm: ${OBJS}
	${CC} ${LFLAGS} ${OBJS} -o qnorm

qnorm.o: qnorm.c
	${CC} ${CFLAGS} Qnorm.c

qfunc.o: qfunc.c 
	${CC} ${CFLAGS} qfunc.c

quicksort.o: quicksort.c 
	${CC} ${CFLAGS} quicksort.c

clean:
	rm -f ${OBJS}
