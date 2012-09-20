CC = cc
CFLAGS	= -Wall -O2 -DHAVE_SINCOS \
	-DPREFIX="./candidates" -DDTAPREFIX="."
LDFLAGS =
LOADLIBES = -lm -lgsl -lgslcblas -lfftw3 -lc

search: search.o JobNAllSky-common.o JobNAllSky.o \
		lvcvirgo.o auxi.o fstat.o gridr.o
sigen:	sigen.o lvcvirgo.o auxi.o
doc:	doxygen ../doc/Doxyfile
clean:
	rm -f *.o

uninstall: clean
	rm -f search sigen
	rm -fr ../doc/html
