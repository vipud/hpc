SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG

include ../config/make.def

OBJS = mg.o ${COMMON}/${RAND}.o ${COMMON}/c_print_results.o \
       ${COMMON}/c_timers.o ${COMMON}/wtime.o \
       ${COMMON}/cl_util.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

mg.o:		mg.c  globals.h npbparams.h
	${CCOMPILE} mg.c

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.h core



