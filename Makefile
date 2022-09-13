
FC = gfortran
FLAGS = -O3  # -fbounds-check
TARGET = lp 

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

LAPACKDIR = /home/gdp/Computazionale/lapack-3.8.0
LAPACK = -L$(LAPACKDIR) -llapack -lrefblas

LIBS = $(LAPACK)

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm *.o *.mod $(TARGET)



interpolators.o: precision.o
functions.o: precision.o interpolators.o
solvers.o: precision.o functions.o
lp.o: precision.o functions.o interpolators.o solvers.o

