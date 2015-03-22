CFLAGS=-std=c99

default: diag

include ${SLEPC_DIR}/conf/slepc_common

diag: diag.o chkopts
	gcc -std=c99  -c qbasic.c vibronic.c
	-${CLINKER} -std=c99 -o diag diag.o qbasic.o vibronic.o ${SLEPC_EPS_LIB} 


main: main.o 
	gcc -std=c99 -c qbasic.c vibronic.c main.c
	gcc -std=c99 -o main main.o qbasic.o vibronic.o
	rm *.o



