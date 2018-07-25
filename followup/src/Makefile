SEARCH_DIR = ../../search/network/src-cpu

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
FFTW_DIR = /opt/fftw/3.3.7-gcc7.2-mvapich2.2/
SINCOS =  YEPPP #YEPPP or SLEEF or GNUSINCOS or NOSINCOS

CC = gcc

# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
CFLAGS = -DPREFIX="./candidates" -DTIMERS -D$(SINCOS)
CFLAGS += -I$(FFTW_DIR)/include -I$(SEARCH_DIR)
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	-O3 -ffast-math -funsafe-loop-optimizations -funroll-loops \
	-march=haswell -mtune=haswell -mavx2 -fopenmp -mfma
#cpu specific optimizations e.g.: 
#sandy-bridge (psk)
	#-march=corei7-avx -mtune=corei7-avx -mavx
#westmere (plgrid)
	#-march=corei7 -mtune=corei7

LDFLAGS = -L$(FFTW_DIR)/lib

LOADLIBES = -Wl,-Bdynamic -lc -lrt -lm -fopenmp -lfftw3 -lgsl -lgslcblas

ifeq ($(strip $(SINCOS)),SLEEF)
SPATH = $(SEARCH_DIR)/lib/sleef-2.80/simd
CFLAGS  += -DENABLE_AVX -I$(SPATH)
LDFLAGS += -L$(SPATH)
LOADLIBES += -lsleef-avx
VECLIB = $(SPATH)/libsleef-avx.so
#sleef-2.80/purec/sleefdp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
VECPATH = $(SEARCH_DIR)/lib/yeppp-1.0.0
CFLAGS  += -I$(VECPATH)/library/headers
LDFLAGS += -L$(VECPATH)/binaries/linux/x86_64
LOADLIBES += -lyeppp
VECLIB = $(VECPATH)/binaries/linux/x86_64/libyeppp.so
endif


all: veclib followup

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	$(CC) -c $(CFLAGS) -Wno-attributes -fPIC $(SPATH)/sleefsimddp.c -o $(SPATH)/sleefsimddp.o
	$(CC) -shared -Wl,-soname,libsleef-avx.so -o $(SPATH)/libsleef-avx.so $(SPATH)/sleefsimddp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
veclib:
	@echo "Using prebuild YEPPP library: $(VECLIB)"
endif


# openmp library added!
followup: followup.o settings.o auxi.o init.o 
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ $(LOADLIBES)

settings.o: settings.c
	$(CC) $(CFLAGS) -c -o $@ $<

init.o: init.c  settings.h struct.h init.h $(SEARCH_DIR)/auxi.h
	$(CC) $(CFLAGS) -c -o $@ $<

auxi.o: $(SEARCH_DIR)/auxi.c $(SEARCH_DIR)/auxi.h
	$(CC) -c -g -O3 -Wall -o $@ $<


# neigh.o: $(SEARCH_DIR)/neigh.c $(SEARCH_DIR)/neigh.h
#	$(CC) $(CFLAGS) -c -o $@ $<



doc:
	doxygen ../doc/Doxyfile

clean:
	rm -f *.o
ifeq ($(strip $(SINCOS)),SLEEF)
	rm -f $(VECLIB)
endif

uninstall: clean
	rm -f search

.PHONY: all veclib
