# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
FFTW_DIR = /opt/fftw/3.3.4/gnu
SINCOS =  YEPPP #YEPPP or SLEEF or GNUSINCOS or NOSINCOS

CC = gcc

# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
CFLAGS = -DPREFIX="./candidates" -DTIMERS -D$(SINCOS)
CFLAGS = -I$(FFTW_DIR)/include 
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	  -O3 -ffast-math -funsafe-loop-optimizations -funroll-loops \
	  -march=corei7-avx -mtune=corei7-avx -mavx
#last line contains cpu specific optimizations e.g.: 
#sandy-bridge (psk)
	#-march=corei7-avx -mtune=corei7-avx -mavx
#westmere (plgrid)
	#-march=corei7 -mtune=corei7

LDFLAGS = -L$(FFTW_DIR)/lib

LOADLIBES = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	    -static -lfftw3 -lgsl -lgslcblas \
	    -Wl,-Bdynamic -lc -lrt -lm

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


all: veclib search

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	$(CC) -c $(CFLAGS) -Wno-attributes -fPIC $(SPATH)/sleefsimddp.c -o $(SPATH)/sleefsimddp.o
	$(CC) -shared -Wl,-soname,libsleef-avx.so -o $(SPATH)/libsleef-avx.so $(SPATH)/sleefsimddp.o
endif
ifeq ($(strip $(SINCOS)),YEPPP)
veclib:
	@echo -e "\nUsing prebuild YEPPP library from: $(VECLIB)\n"
endif


search: main.o jobcore.o timer.o settings.o auxi.o init.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LOADLIBES)

jobcore.o: settings.h timer.h auxi.h

main.o: settings.h struct.h init.h auxi.h

init.o: settings.h struct.h init.h auxi.h

settings.o: settings.h struct.h

auxi.o: auxi.c auxi.h
	$(CC) -c -O3 -Wall -o $@ $< 

clean:
	rm -f *.o
ifeq ($(strip $(SINCOS)),SLEEF)
	rm -f $(VECLIB)
endif

uninstall: clean
	rm -f search

.PHONY: all veclib clean uninstall

#----------------------------------------------

#sigen: sigen.o settings.o auxi.o
#	$(CC) $(LDFLAGS) -o $@ $^ $(LOADLIBES)
#
#gridopt.o: gridopt.c 
#	$(CC) -c -o $@ $< -O3 -Wall
#
#lineph.o: lineph.c 
#	$(CC) -c -o $@ $< -O3 -Wall
#
#gridopt: gridopt.o auxi.o lineph.o settings.o
#	$(CC) $(LDFLAGS) -o $@ $^ -Wl,-Bdynamic -lc -lrt -lgsl -lgslcblas -lm 
