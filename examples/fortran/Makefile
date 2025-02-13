FC = gfortran
FFLAGS = -Wall -Wextra -O0 -g -ffree-line-length-none -Wno-unused-dummy-argument
PINEAPPL_LIBS != pkg-config pineappl_capi --libs
LHAPDF_LIBS != pkg-config lhapdf --libs

PROGRAMS = \
	dyaa-deprecated \
	test \
	test-deprecated \
	lhapdf_example \
	lhapdf_example-deprecated

all: $(PROGRAMS)

test-examples: $(PROGRAMS)
	set -e && for i in $(PROGRAMS); do ./$${i}; done > test-output && diff -u output test-output && rm test-output

pineappl.o: pineappl.f90
	$(FC) $(FFLAGS) -c $< -o $@

dyaa-deprecated: dyaa-deprecated.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(PINEAPPL_LIBS) -o $@

test: test.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(PINEAPPL_LIBS) -o $@

test-deprecated: test-deprecated.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(PINEAPPL_LIBS) -o $@

lhapdf_example: lhapdf_example.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o $@

lhapdf_example-deprecated: lhapdf_example-deprecated.f90 pineappl.o
	$(FC) $(FFLAGS) $< pineappl.o $(LHAPDF_LIBS) $(PINEAPPL_LIBS) -o $@

PHONY: clean

clean:
	rm -f *.o *.mod $(PROGRAMS)
