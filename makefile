default: diag

include ${SLEPC_DIR}/conf/slepc_common

diag: diag.o chkopts
	gcc -c qbasic.c vibronic.c
	-${CLINKER} -o diag diag.o qbasic.o vibronic.o ${SLEPC_EPS_LIB}

