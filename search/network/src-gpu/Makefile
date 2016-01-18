TARGET = gwsearch-gpu

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
FFTW_DIR = /opt/fftw/3.3.4/gnu
SINCOS =  GNUSINCOS #or NOSINCOS

CC = gcc
NVCC = gcc -x c #nvcc

CFLAGS = -DPREFIX="./candidates" -DTIMERS=1 -D$(SINCOS) -DVERBOSE
NCFLAGS := $(CFLAGS)
CFLAGS += -I$(FFTW_DIR)/include 
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
          -O3 -ffast-math -funsafe-loop-optimizations -funroll-loops \
          -march=corei7-avx -mtune=corei7-avx -mavx
#last line contains cpu specific optimizations e.g.: 
#sandy-bridge (psk)
        #-march=corei7-avx -mtune=corei7-avx -mavx
#westmere (plgrid)
        #-march=corei7 -mtune=corei7
LDFLAGS = -L$(FFTW_DIR)/lib
LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
         -static -lfftw3 -lgsl -lgslcblas \
         -Wl,-Bdynamic -lc -lrt -lm

NCFLAGS += -I$(FFTW_DIR)/include
NCFLAGS += -O3 -g #--ptxas-options=-v -G -use_fast_math -arch=sm_35
NCFLAGS += -I/usr/local/cuda/include/


#MACROS = -DPREFIX="./candidates" -DDTAPREFIX="." -DTIMERS=1 -DGNUSINCOS
#CFLAGS = -O3 -g
#NCFLAGS = --ptxas-options=-v -G -use_fast_math -arch=sm_35 
#
#INC = -I/usr/local/cuda/include/
#
#LDFLAGS = -L/usr/local/cuda/lib64
#LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
#	 -static -lgsl -lgslcblas \
#	 -Wl,-Bdynamic -lcufft -lcuda -lcudart -lcusparse -lcublas -lc -lrt -lm
#

OBJS = main.o jobcore.o timer.o settings.o auxi.o init.o #spline_z.o kernels.o

${TARGET}: ${OBJS}
	$(CC) ${LDFLAGS} -o $@ $+ $(LDLIBS)

%.o: %.c
	$(CC) -c -Wall ${CFLAGS} -o $@ $+

%.o: %.cu
	$(NVCC) -c $(NCFLAGS) -o $@ $+

clean:
	rm -f ${TARGET} ${OBJS}

.PHONY: $(TARGET) clean
