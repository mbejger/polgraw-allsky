TARGET=coincidences
SEARCH_DIR=../../search/network/src-cpu

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
FFTW_DIR = /opt/fftw/3.3.5/gnu-omp-devtoolset3

CC = gcc -std=gnu11 -static-libgcc
# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
#CFLAGS = -DPREFIX="./candidates" -DTIMERS
CFLAGS = -I$(SEARCH_DIR) -I$(FFTW_DIR)/include
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	  -O3 -funsafe-loop-optimizations -funroll-loops \
		-march=broadwell -mtune=broadwell -mavx2
#	  -march=sandybridge -mtune=sandybridge -mavx

#cpu specific optimizations e.g.: 
#haswell
#	  -march=haswell -mtune=haswell -mavx2 -mfma
#sandy-bridge (psk)
	#-march=sandybridge -mtune=sandybridge -mavx
#westmere (plgrid)
	#-march=corei7 -mtune=corei7

LDFLAGS = -L$(FFTW_DIR)/lib

LOADLIBES = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	    -static -lfftw3 -lgsl -lgslcblas \
	    -Wl,-Bdynamic -lc -lrt -lm


all: $(TARGET)

$(TARGET): coincidences.o settings.o auxi.o init.o
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) -o $@

coincidences.o: coincidences.c

settings.o: $(SEARCH_DIR)/settings.c
	$(CC) $(CFLAGS) -c -o $@ $<

init.o: $(SEARCH_DIR)/init.c  $(SEARCH_DIR)/settings.h $(SEARCH_DIR)/struct.h $(SEARCH_DIR)/init.h $(SEARCH_DIR)/auxi.h
	$(CC) $(CFLAGS) -c -o $@ $<

auxi.o: $(SEARCH_DIR)/auxi.c $(SEARCH_DIR)/auxi.h
	$(CC) $(CFLAGS) -c -o $@ $<
#$(CC) -c -g -O3 -Wall -o $@ $<

fap: fap.o settings.o auxi.o
	$(CC) -o $@ $^ -lc -lrt -lm -lgsl -lgslcblas -lfftw3

fap-many: fap-many.o settings.o auxi.o
	$(CC) -o $@ $^ -lc -lrt -lm -lgsl -lgslcblas -lfftw3


doc:
	doxygen ../doc/Doxyfile

clean:
	rm -f *.o

uninstall: clean
	rm -f coincidences

.PHONY: all doc clean uninstall
