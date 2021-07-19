SHELL=/bin/sh
OBJECTS=main.o riemann.o matrix.o silo.o mesh.o shape.o quadrature.o polynomial.o
MYFLAGS=-g -O2 -Wall -v -std=c++11 -pthread
MYCMP=g++
MYINCS=-I/usr/include/hdf5/serial
MYLIBS=
EXEC=rholo
BLAS=-L/usr/lib/x86_64-linux-gnu -lblas
LAPACK=-L/usr/lib/x86_64-linux-gnu -llapack
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
SILO=/usr/lib/libsiloh5.a

${EXEC}: ${OBJECTS}
	${MYCMP} ${MYFLAGS} ${MYINCS} -o $@ ${OBJECTS} ${SILO} ${HDF5} -lz -ldl -lsz ${MYLIBS} ${LAPACK} ${BLAS}

nolink: ${OBJECTS}
	${MYCMP} ${MYFLAGS} -c ${OBJECTS}
	-echo ${OBJECTS}

clean:
	-rm -f *.o ${EXEC}

veryclean: clean
	-rm -f ../bin/${EXEC}

.cpp.o:
	${MYCMP} ${MYFLAGS} ${MYINCS} -c $<
