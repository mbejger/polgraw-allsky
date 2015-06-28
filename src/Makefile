CC = gcc
CFLAGS	= -Wall -O3 -DPREFIX="./candidates" -DDTAPREFIX="." -DTIMERS=1
LDFLAGS =
LIBS = -lm -lgsl -lgslcblas -lfftw3 -lc -lrt

search: main.o\
  jobcore.o timer.o settings.o auxi.o init.o 
	$(CC) -gnatv -o $@ $^ $(LIBS)

sigen: sigen.o\
	settings.o auxi.o
	$(CC) -gnatv -o $@ $^ $(LIBS)

gridopt: gridopt.o\
	auxi.o lineph.o settings.o
	$(CC) -gnatv -o $@ $^ $(LIBS)

clean:
	rm -f *.o

uninstall: clean
	rm -f search sigen gridopt
