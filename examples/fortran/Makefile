FC = gfortran
FFLAGS = -Wall -Wextra -O3
PINEAPPL_LIBS != pkg-config pineappl_capi --libs

all: pineappl.f90 dyaa.f90
	$(FC) $(FFLAGS) pineappl.f90 dyaa.f90 $(PINEAPPL_LIBS) -o dyaa
