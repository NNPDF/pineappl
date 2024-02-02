FC = gfortran
FFLAGS = -Wall -Wextra -O0 -g
PINEAPPL_LIBS != pkg-config pineappl_capi --libs
LHAPDF_LIBS != pkg-config lhapdf --libs

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: pineappl.o dyaa.o test.o
	$(FC) $(FFLAGS) dyaa.o pineappl.o $(PINEAPPL_LIBS) -o dyaa
	$(FC) $(FFLAGS) test.o pineappl.o $(PINEAPPL_LIBS) -o test

lhapdf_example: pineappl.o lhapdf_example.o
	$(FC) $(FFLAGS) lhapdf_example.o pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o lhapdf_example

clean:
	rm -f *.o *.mod dyaa test lhapdf_example
