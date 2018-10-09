SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG

include ../config/make.def

include ../sys/make.common

OBJS = cg.o ${COMMON}/${RAND}.o ${COMMON}/c_print_results.o \
       ${COMMON}/c_timers.o ${COMMON}/wtime.o \
       ${COMMON}/cl_util.o

${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

.c.o:
	${CCOMPILE} $<

cg.o: cg.c  mpinpb.h npbparams.h 

clean:
	- rm -f *.o *~ mputil* *.bin
	- rm -f cg npbparams.h core
	- if [ -d rii_files ]; then rm -r rii_files; fi
