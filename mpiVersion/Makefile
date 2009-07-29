CC = mpicc
LFLAGS = -lm

OBJS = Qnorm.o Qfunc.o quicksort.o

Qnorm: ${OBJS}
	${CC} ${LFLAGS} ${OBJS} -o Qnorm

qnorm.o: qnorm.c
	${CC} ${LFLAGS} -c Qnorm.c

qFunc.o: Qfunc.c 
	${CC} ${LFLAGS} -c qFunc.c

quicksort.o: quicksort.c 
	${CC} ${LFLAGS} -c quicksort.c

clean:
	rm -f ${OBJS}
