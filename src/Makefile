MACROS = -DPREFIX="./candidates" -DDTAPREFIX="." -DTIMERS=1
CFLAGS = -O3 -g -I/usr/local/cuda/include/
LDFLAGS = -L/usr/local/cuda/lib -L/usr/local/cuda/lib64
LIBS = -lm -lgsl -lgslcblas -lcufft -lcuda -lcudart -lcusparse

TARGET = gwsearch-gpu
OBJS = main.o jobcore.o auxi.o timer.o settings.o init.o spline_z.o kernels.o

${TARGET}: ${OBJS}
	g++ ${LDFLAGS} -o $@ $+ $(LIBS)

%.o: %.c
	g++ -c -Wall -pedantic -std=c++0x ${CFLAGS} ${MACROS} -o $@ $+

%.o: %.cu
	nvcc --ptxas-options=-v -c -G --use_fast_math -arch=sm_35 ${CFLAGS} ${MACROS} -o $@ $+

clean:
	rm -f ${TARGET} ${OBJS}
