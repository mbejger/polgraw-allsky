.PHONY: all veclib
# adjust FFTW_DIR and SINCOS
FFTW_DIR = /opt/fftw/3.3.4/gnu
SINCOS =  YEPPP #YEPPP or SLEEF or GNUSINCOS or NOSINCOS

CC = gcc -g -fno-omit-frame-pointer
CFLAGS = -Wall -Wno-unused -O3 -DPREFIX="./candidates" -DTIMERS=1 \
	-march=native -mtune=native \
	-ffast-math -funroll-loops -funsafe-loop-optimizations \
        -I$(FFTW_DIR)/include -D$(SINCOS)
LDFLAGS = -L$(FFTW_DIR)/lib
LOADLIBES = -lm -lgsl -lgslcblas -lfftw3 -lc

ifeq ($(strip $(SINCOS)),SLEEF)
SPATH = lib/sleef-2.80/simd
CFLAGS  += -DENABLE_AVX -I$(SPATH)
LDFLAGS += -L$(SPATH)
LOADLIBES += -lsleef-avx
VECLIB = $(SPATH)/libsleef-avx.so
#sleef-2.80/purec/sleefdp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
CFLAGS  += -mavx -Ilib/yeppp-1.0.0/library/headers
LDFLAGS += -Llib/yeppp-1.0.0/binaries/linux/x86_64
LOADLIBES +=  -lyeppp
VECLIB = lib/yeppp-1.0.0/binaries/linux/x86_64/libyeppp.so
endif


all: veclib search

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	$(CC) -c $(CFLAGS) -Wno-attributes -fPIC $(SPATH)/sleefsimddp.c -o $(SPATH)/sleefsimddp.o
	$(CC) -shared -Wl,-soname,libsleef-avx.so -o $(SPATH)/libsleef-avx.so $(SPATH)/sleefsimddp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
veclib:
	echo "Using prebuild YEPPP library: $(VECLIB)"
endif

search: search.o JobNAllSky-common.o JobNAllSky.o \
	settings.o auxi.o fstat.o gridr.o

bench: bench.o settings.o

doc:
	doxygen ../doc/Doxyfile

gridopt: gridopt.o\
	auxi.o lineph.o settings.o
	$(CC) -gnatv -o $@ $^ $(LIBS)

clean:
	rm -f *.o
ifeq ($(strip $(SINCOS)),SLEEF)
	rm -f $(VECLIB)
endif

uninstall: clean
	rm -f search bench 
	rm -fr ../doc/html
