CC = g++
NVCC = nvcc

MACROS = -DPREFIX="./candidates" -DDTAPREFIX="." -DTIMERS=1
CFLAGS = -O3 -g
NCFLAGS = --ptxas-options=-v -G -use_fast_math -arch=sm_35 

INC = -I/usr/local/cuda/include/

LDFLAGS = -L/usr/local/cuda/lib -L/usr/local/cuda/lib64
LIBS = -static -lgsl -lgslcblas -Wl,-Bdynamic -lcufft -lcuda -lcudart -lcusparse -lcublas -lc -lrt -lm

TARGET = gwsearch-gpu
OBJS = main.o jobcore.o auxi.o timer.o settings.o init.o spline_z.o kernels.o

${TARGET}: ${OBJS}
	$(CC) ${LDFLAGS} -o $@ $+ $(LIBS)

%.o: %.c
	$(CC) -c -Wall -pedantic -std=c++0x ${CFLAGS} ${MACROS} $(INC) -o $@ $+

%.o: %.cu
	$(NVCC) -c ${CFLAGS} $(NCFLAGS) ${MACROS} $(INC) -o $@ $+

clean:
	rm -f ${TARGET} ${OBJS}

.PHONY: $(TARGET) clean
