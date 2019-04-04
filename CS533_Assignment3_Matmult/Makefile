# GCC Compiler Flags
CC = gcc 
OPT = -O3 -mfpmath=sse -msse4.2 -g -funroll-all-loops -ftree-vectorize -ftree-slp-vectorize -fopenmp-simd -ftree-vectorizer-verbose=2
CFLAGS += -Wall -std=gnu11 -DGETTIMEOFDAY $(OPT) -I/opt/ohpc/pub/libs/gnu/gsl/2.4/include
LDLIBS = -L/opt/ohpc/pub/libs/gnu/gsl/2.4/lib -lgslcblas

# Intel compiler flags
#CC = icc
#OPT = -O3 -qopenmp
#CFLAGS += -std=c99 -Wall -DGETTIMEOFDAY $(OPT)
#LDLIBS =  -mkl


types = naive basic row col blocked rb copy blas genvect autovect
#vect vect_copy vect_copy_2
targets = $(patsubst %,benchmark_%,$(types))
objects = benchmark.o $(patsubst %,dgemm_%.o,$(types))

.PHONY : default
default : all

.PHONY : all
all : $(targets) 

objects : $(objects)

benchmark_% : benchmark.o dgemm_%.o 
	$(CC) $(OPT) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects) *.stdout

%.o : %.c
	$(CC) -c $(CFLAGS) $<
