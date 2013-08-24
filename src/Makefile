OPT=-O3 -march=native #-pg
CC=gcc $(OPT)
#NCC=nvcc -keep --ptxas-options=-v -arch sm_20 -O3 -Xptxas -dlcm=cg
NCC=nvcc -arch sm_20 -O3

CPP_DEFS=-DTIMERS=1
CPP_DEFS_CU=-DCUDA_DEV=0 -DBLOCK_DIM=128 $(CPP_DEFS)
CFLAGS=-I.

# librt is needed by clock_gettime only
LIBS=-lm -lrt

DEPS = 
OBJ = search.o settings.o gridr.o auxi.o

#this in needed to avoid shitstorms when test.cu.cpp is left after previous make
%.cu:

%.o: %.cu $(DEPS)
	$(NCC) $(CPP_DEFS_CU) $(CFLAGS) -o $@ -c $< 

%.o: %.c $(DEPS)
	$(CC) $(CPP_DEFS) $(CFLAGS) -o $@ -c $< 

test: search.o settings.o gridr.o auxi.o
	$(NCC) -o $@ $^ $(LIBS)

all: search

clean:
	rm -f $(OBJ)
	$(NCC) --clean-targets -c search.cu -o search.o
