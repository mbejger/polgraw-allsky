CC = gcc
CFLAGS	= -Wall -O3 -DHAVE_SINCOS 
LDFLAGS =
LOADLIBES = -lm -lgsl -lgslcblas -lfftw3 -lc

search: search.o JobNAllSky-common.o JobNAllSky.o \
	settings.o auxi.o fstat.o gridr.o
bench: bench.o settings.o
doc:
	doxygen ../doc/Doxyfile
clean:
	rm -f *.o

uninstall: clean
	rm -f search bench 
	rm -fr ../doc/html
