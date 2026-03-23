FC = gfortran
FFLAGS = -Wall -Wextra -O0 -g -ffree-line-length-none
PINEAPPL_LIBS != pkg-config pineappl_capi --libs
LHAPDF_LIBS != pkg-config lhapdf --libs

PROGRAMS = \
	dyaa \
	test \
	lhapdf_example

all: $(PROGRAMS)

test-examples: $(PROGRAMS)
	set -e && for i in $(PROGRAMS); do ./$${i} > output; diff -u $${i}.output output; done; rm output

pineappl.o: pineappl.f90
	$(FC) $(FFLAGS) -c $< -o $@

dyaa: dyaa.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(PINEAPPL_LIBS) -o $@

test: test.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(PINEAPPL_LIBS) -o $@

lhapdf_example: lhapdf_example.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o $@

.PHONY: clean

clean:
	rm -f *.o *.mod $(PROGRAMS)
