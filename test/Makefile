
FLAGS 		= 

HEALPIX_LIBS = -L/Users/chervias/Software/miniconda3/envs/cmb/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
#HEALPIX_LIBS = -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
HEALPIX_INCLUDE = -I/Users/chervias/Software/miniconda3/envs/cmb/include/healpix_cxx

OBJS = test.o

all: test

test: ${OBJS}
	g++ ${OBJS} -o test ${FLAGS} -lm $(HEALPIX_LIBS)

test.o: test.cpp
	g++ -o ${OBJS} ${HEALPIX_INCLUDE} -c test.cpp ${FLAGS}

clean:
	rm -rf test.o test
