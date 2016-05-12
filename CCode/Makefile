COMP = gcc
LIB = -lm -lgsl -lgslcblas
OPT = -Wall -O3 -DHAVE_INLINE
OBJECTS = main.o partition.o SearchAlgs.o

FindGroups: ${OBJECTS}
	${COMP} ${OBJECTS} -o FindGroups ${OPT} ${LIB}

main.o: main.c Common.h SearchAlgs.h partition.h
	${COMP} main.c -c ${OPT}

partition.o: partition.c Common.h
	${COMP} partition.c -c ${OPT}

SearchAlgs.o: SearchAlgs.c partition.h Common.h
	${COMP} SearchAlgs.c -c ${OPT}

clean: 
	rm -rf *.o *~ *.log FindGroups
