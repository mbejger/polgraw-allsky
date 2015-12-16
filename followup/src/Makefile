SEARCH_DIR = ../../search/network/src-cpu

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
FFTW_DIR = /opt/fftw/3.3.4/gnu
SINCOS =  YEPPP #YEPPP or SLEEF or GNUSINCOS or NOSINCOS

CC = gcc

# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
CFLAGS = -DPREFIX="./candidates" -DTIMERS -D$(SINCOS)
CFLAGS += -I$(FFTW_DIR)/include -I$(SEARCH_DIR)
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	-O3 -ffast-math -funsafe-loop-optimizations -funroll-loops \
	-march=corei7-avx -mtune=corei7-avx -mavx -fopenmp 
#cpu specific optimizations e.g.: 
#sandy-bridge (psk)
	#-march=corei7-avx -mtune=corei7-avx -mavx
#westmere (plgrid)
	#-march=corei7 -mtune=corei7

LDFLAGS = -L$(FFTW_DIR)/lib

LOADLIBES = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	    -static -lfftw3 -lgsl -lgslcblas \
	    -Wl,-Bdynamic -lc -lrt -lm -fopenmp 

ifeq ($(strip $(SINCOS)),SLEEF)
SPATH = lib/sleef-2.80/simd
CFLAGS  += -DENABLE_AVX -I$(SPATH)
LDFLAGS += -L$(SPATH)
LOADLIBES += -lsleef-avx
VECLIB = $(SPATH)/libsleef-avx.so
#sleef-2.80/purec/sleefdp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
CFLAGS  += -Ilib/yeppp-1.0.0/library/headers
LDFLAGS += -Llib/yeppp-1.0.0/binaries/linux/x86_64
LOADLIBES += -lyeppp
VECLIB = lib/yeppp-1.0.0/binaries/linux/x86_64/libyeppp.so
endif


all: veclib followup

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	$(CC) -c $(CFLAGS) -Wno-attributes -fPIC $(SPATH)/sleefsimddp.c -o $(SPATH)/sleefsimddp.o
	$(CC) -shared -Wl,-soname,libsleef-avx.so -o $(SPATH)/libsleef-avx.so $(SPATH)/sleefsimddp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
veclib:
	echo "Using prebuild YEPPP library: $(VECLIB)"
endif


# open mp library added!
followup: followup.o settings.o auxi.o init.o 
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ $(LOADLIBES)

settings.o: $(SEARCH_DIR)/settings.c
	$(CC) $(CFLAGS) -c -o $@ $<

init.o: $(SEARCH_DIR)/init.c  $(SEARCH_DIR)/settings.h $(SEARCH_DIR)/struct.h $(SEARCH_DIR)/init.h $(SEARCH_DIR)/auxi.h
	$(CC) $(CFLAGS) -c -o $@ $<

auxi.o: $(SEARCH_DIR)/auxi.c $(SEARCH_DIR)/auxi.h
	$(CC) -c -g -O3 -Wall -o $@ $<


doc:
	doxygen ../doc/Doxyfile

#auxi.o: auxi.c 
#	$(CC) -c -o $@ $< -O3 -Wall

clean:
	rm -f *.o
ifeq ($(strip $(SINCOS)),SLEEF)
	rm -f $(VECLIB)
endif

uninstall: clean
	rm -f search

.PHONY: all veclib
