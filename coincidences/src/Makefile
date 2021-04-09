#TARGET=coincidences
TARGET=coincidences-haswell-dev7
SEARCH_DIR=../../search/network/src-openmp

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
FFTW_DIR = /opt/fftw/3.3.7-gcc7.2-mvapich2.2

CC = gcc -std=gnu11 -static-libgcc
# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
#CFLAGS = -DPREFIX="./candidates" -DTIMERS
CFLAGS = -I$(SEARCH_DIR) -I$(FFTW_DIR)/include $(shell gsl-config --cflags)
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	  -O3 -funsafe-loop-optimizations -funroll-loops \
	  -march=haswell -mtune=haswell -mavx2 -mfma

#cpu specific optimizations e.g.: 
#haswell
	#-march=haswell -mtune=haswell -mavx2 -mfma
#sandy-bridge (psk)
	#-march=sandybridge -mtune=sandybridge -mavx
#westmere (plgrid)
	#-march=corei7 -mtune=corei7

LDFLAGS=-L$(FFTW_DIR)/lib -L$(shell gsl-config --prefix)/lib

LOADLIBES=-Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
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

fap: fap.o settings.o auxi.o
	$(CC) -o $@ $^ -lc -lrt -lm -lgsl -lgslcblas -lfftw3

fap_new_coi_header: fap_new_coi_header.o settings.o auxi.o
	$(CC) -o $@ $^ -lc -lrt -lm -lgsl -lgslcblas -lfftw3

fap_new_coi_header_hyper_improved: fap_new_coi_header_hyper_improved.o settings.o auxi.o
	$(CC) -o $@ $^ -lc -lrt -lm -lgsl -lgslcblas -lfftw3

fap2: fap2.o settings.o auxi.o
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) -o $@

doc:
	doxygen ../doc/Doxyfile

clean:
	rm -f *.o

uninstall: clean
	rm -f coincidences

.PHONY: all doc clean uninstall
