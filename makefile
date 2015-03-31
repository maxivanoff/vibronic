CFLAGS=-std=c99

default:diag

include ${SLEPC_DIR}/conf/slepc_common

diag: diag.o chkopts
	gcc -std=c99  -c multi.c
	-${CLINKER} -std=c99 -o diag diag.o multi.o ${SLEPC_EPS_LIB} 


main: main.o chkopts
	gcc -std=c99 -c  multi.c main.c
	-${CLINKER} -std=c99 -o main main.o  multi.o
	rm *.o



