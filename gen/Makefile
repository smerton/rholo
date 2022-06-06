SHELL=/bin/sh

OBJECTS=input.o output.o main.o
MYFLAGS=-g -O2 -Wall -v -std=c++11 -pthread
MYCMP=g++
MYINCS=
MYLIBS=
EXEC=gen

${EXEC}: ${OBJECTS}
	${MYCMP} ${MYFLAGS} ${MYINCS} -o $@ ${OBJECTS} ${MYLIBS}
	-cp ${EXEC} ../bin

nolink: ${OBJECTS}
	${MYCMP} ${MYFLAGS} -c ${OBJECTS}
	-echo ${OBJECTS}

clean:
	-rm -f *.o ${EXEC}

veryclean: clean
	-rm -f ../bin/${EXEC}

.cpp.o:
	${MYCMP} ${MYFLAGS} ${MYINCS} -c $<
