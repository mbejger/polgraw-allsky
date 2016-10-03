TARGET = gwsearch-gpu

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
FFTW_DIR = /opt/fftw/3.3.4/gnu
SINCOS =  GNUSINCOS #or NOSINCOS

CC = g++
NVCC = nvcc -rdc=true

CFLAGS = -DPREFIX="./candidates" -DTIMERS=3 -D$(SINCOS) -DVERBOSE -DCUDA_DEV=0
CFLAGS += -I/usr/local/cuda/include/
NCFLAGS := $(CFLAGS)
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
          -O3 -ffast-math -funsafe-loop-optimizations -funroll-loops \
          -march=native -mtune=native -mavx
#last line contains cpu specific optimizations e.g.: 
#sandy-bridge (psk)
        #-march=corei7-avx -mtune=corei7-avx -mavx
#westmere (plgrid)
        #-march=corei7 -mtune=corei7
#LDFLAGS = -L$(FFTW_DIR)/lib
#LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
#         -static -lfftw3 -lgsl -lgslcblas \
#         -Wl,-Bdynamic -lc -lrt -lm

#NCFLAGS += -I$(FFTW_DIR)/include 
NCFLAGS += -DBLOCK_SIZE=256 
NCFLAGS += -O3 -g --ptxas-options=-v -G -use_fast_math -arch=sm_35


LDFLAGS = -L/usr/local/cuda/lib64
#LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
#	 -static -lgsl -lgslcblas \
#	 -Wl,-Bdynamic -lcufft -lcuda -lcudart -lcusparse -lcublas -lc -lrt -lm

LDLIBS = -Xlinker --dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	 -Xlinker -static -lgsl -lgslcblas \
	 -Xlinker -Bdynamic -lcufft -lcuda -lcudart -lcusparse -lcublas -lc -lrt -lm

OBJS = main.o jobcore.o timer.o settings.o auxi.o init.o kernels.o spline_z.o

${TARGET}: ${OBJS}
	$(NVCC) ${NCFLAGS} ${LDFLAGS} -o $@ $+ $(LDLIBS)

%.o: %.c
	$(CC) -c -Wall ${CFLAGS} -o $@ $+

%.o: %.cu
	$(NVCC) -c $(NCFLAGS) -o $@ $+

clean: 
	$(NVCC) --clean-targets -c jobcore.cu -o jobcore.o
	$(NVCC) --clean-targets -c init.cu -o init.o
	$(NVCC) --clean-targets -c kernels.cu -o kernels.o
	$(NVCC) --clean-targets -c spline_z.cu -o spline_z.o
	rm -f ${TARGET} ${OBJS}


.PHONY: $(TARGET) clean
