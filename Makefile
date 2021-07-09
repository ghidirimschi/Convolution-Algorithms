CC = gcc
CFLAGS  = -Wall -O2

OBJECTS = timer.o LibList.o main.o direct.o fftconv.o karatsuba.o dupl.o sietse.o lrsno.o lrsnosort.c oamfft.o derivative.o nttbarret.o nttsoftware.o nttcrt.o nopadfft.o nttcrtbarrett.o

all: ${OBJECTS}
	${CC} ${CLAGS} -o main ${OBJECTS} -lm
tar:
	tar cvf convolution.tar *.[ch] Makefile
clean:
	rm -f *~
	rm -f *.o
