CC = cc
#LAL_FLAGS = -DUSE_LAL -I/usr/share/lalpulsar/include -I/usr/include/lal/include
LAL_FLAGS = 
#LAL_LIBS = -llalpulsar -L/usr/share/lalpulsar/lib
LAL_LIBS =
INIPARSER_DIR = iniparser
INIPARSER_FLAGS = -I$(INIPARSER_DIR)/src
INIPARSER_LIBS = -Liniparser -liniparser
CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -O2 -DDEBUG
# debug flags below
#CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -DDEBUG -g
LOADLIBES = -lm -lfftw3 -lgsl -lgslcblas $(LAL_LIBS) $(INIPARSER_LIBS)

.PHONY:	iniparser

all:	iniparser gen2day genseg

iniparser:	
	$(MAKE) -C $(INIPARSER_DIR)

gen2day: ExtractBand.o DeleteOutliers.o get_barycenter.o sid.o
	$(CC) -o gen2day  $^ $(LOADLIBES)

genseg:	genseg.o GrubbsOutliers.o get_barycenter.o sid.o
	$(CC) -o genseg  $^ $(LOADLIBES)

laltest: laltest.o get_barycenter.o sid.o
	$(CC) -o laltest $^ $(LAL_LIBS)

clean:
	rm -f *.o gen2day laltest
	$(MAKE) -C $(INIPARSER_DIR) veryclean
