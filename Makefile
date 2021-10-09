SHELL=/bin/sh
#OBJECTS=main.o riemann.o matrix.o silo.o mesh.o shape.o timer.o quadrature.o polynomial.o
OBJECTS=main.o silo.o eos.o riemann.o matrix.o mesh.o shape.o timer.o quadrature.o polynomial.o
MYFLAGS=-g -O2 -Wall -v -std=c++11 -pthread
MYCMP=g++
#MYINCS=-I/usr/include/hdf5/serial -I../../typhonio/master/build/gcc/include
MYINCS=-I/usr/include/hdf5/serial
MYLIBS=
EXEC=rholo
BLAS=-L/usr/lib64 -lblas
LAPACK=-L/usr/lib64 -llapack
#TYPHONIO=-L../../typhonio/master/build/gcc/lib -ltyphonio
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -lz -ldl -lsz
SILO=/usr/lib/libsiloh5.a
ZLIBS= -lz -ldl -lsz

${EXEC}: ${OBJECTS}
	${MYCMP} ${MYFLAGS} ${MYINCS} -o $@ ${OBJECTS} ${SILO} ${TYPHONIO} ${HDF5} ${ZLIBS} ${MYLIBS} ${LAPACK} ${BLAS}

nolink: ${OBJECTS}
	${MYCMP} ${MYFLAGS} -c ${OBJECTS}
	-echo ${OBJECTS}

clean:
	-rm -f *.o ${EXEC}

veryclean: clean
	-rm -f ../bin/${EXEC}

.cpp.o:
	${MYCMP} ${MYFLAGS} ${MYINCS} -c $<
