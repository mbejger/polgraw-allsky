CC = cc
CFLAGS	= -Wall -O3 -DPREFIX="./candidates" -DDTAPREFIX="." -DTIMERS=1
LDFLAGS =
LIBS = -lm -lgsl -lgslcblas -lfftw3 -lc -lrt

search: main.o jobcore.o timer.o\
	settings.o auxi.o init.o 
	$(CC) -gnatv -o $@ $^ $(LIBS)

clean:
	rm -f *.o

uninstall: clean
	rm -f search
