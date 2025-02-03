FC = gfortran
FFLAGS = -Wall -Wextra -O0 -g
PINEAPPL_LIBS != pkg-config pineappl_capi --libs
LHAPDF_LIBS != pkg-config lhapdf --libs

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: pineappl.o dyaa-deprecated.o test.o test-deprecated.o
	$(FC) $(FFLAGS) dyaa-deprecated.o pineappl.o $(PINEAPPL_LIBS) -o dyaa-deprecated
	$(FC) $(FFLAGS) test.o pineappl.o $(PINEAPPL_LIBS) -o test
	$(FC) $(FFLAGS) test-deprecated.o pineappl.o $(PINEAPPL_LIBS) -o test-deprecated

lhapdf_examples: pineappl.o lhapdf_example.o lhapdf_example-deprecated.o
	$(FC) $(FFLAGS) lhapdf_example.o pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o lhapdf_example
	$(FC) $(FFLAGS) lhapdf_example-deprecated.o pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o lhapdf_example-deprecated

clean:
	rm -f *.o *.mod dyaa-deprecated test test-deprecated lhapdf_example lhapdf_example-deprecated
