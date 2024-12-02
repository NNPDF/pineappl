FC = gfortran
FFLAGS = -Wall -Wextra -O0 -g
PINEAPPL_LIBS != pkg-config pineappl_capi --libs
LHAPDF_LIBS != pkg-config lhapdf --libs

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: pineappl.o dyaa.o test.o test_v1.o
	$(FC) $(FFLAGS) dyaa.o pineappl.o $(PINEAPPL_LIBS) -o dyaa
	$(FC) $(FFLAGS) test.o pineappl.o $(PINEAPPL_LIBS) -o test
	$(FC) $(FFLAGS) test_v1.o pineappl.o $(PINEAPPL_LIBS) -o test_v1

lhapdf_examples: pineappl.o lhapdf_example.o lhapdf_example_v1.o
	$(FC) $(FFLAGS) lhapdf_example.o pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o lhapdf_example
	$(FC) $(FFLAGS) lhapdf_example_v1.o pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o lhapdf_example_v1

clean:
	rm -f *.o *.mod dyaa test test_v1 lhapdf_example lhapdf_example_v1
