
FLAGS 		= 

HEALPIX_LIBS = -L/home/chervias/Software/anaconda3/envs/healpy/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
#HEALPIX_LIBS = -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
HEALPIX_INCLUDE = -I/home/chervias/Software/anaconda3/envs/healpy/include/healpix_cxx

OBJS = test.o

all: test

test: ${OBJS}
	g++ ${OBJS} -o test ${FLAGS} -lm $(HEALPIX_LIBS)

test.o: test.cpp
	g++ -o ${OBJS} ${HEALPIX_INCLUDE} -c test.cpp ${FLAGS}

clean:
	rm -rf test.o test
