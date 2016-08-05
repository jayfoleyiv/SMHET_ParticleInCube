LIB        = -L. -lm -lfftw3 -lblas -Lg2c
INCLUDE    = -I.
CFLAGS     = -O2
EXEC       = PIB.x
CXX        = g++

${EXEC}: Particle_in_Cube.c blas.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} Particle_in_Cube.c blas.o -o ${EXEC}

blas.o: blas.c blas.h
	${CXX} ${LIB} -c blas.c ${CFLAGS}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

