CXX = clang++
FC = gfortran

# Ask gsl-config for flags
GSL_CFLAGS  = $(shell gsl-config --cflags)
GSL_LIBS    = $(shell gsl-config --libs)

CXXFLAGS = -Wall -Wfatal-errors -std=c++17 -O3 $(GSL_CFLAGS)
FFLAGS = -O3
LIBS = $(GSL_LIBS) -lstdc++ -lgfortran

# Use rpath so the binary knows where conda’s libgsl is
RPATH = -Wl,-rpath,$(CONDA_PREFIX)/lib

all: GWB_SMBHB

GWB_SMBHB: GWB_SMBHB.o bessel/CBESSEL.o  bessel/fortran_bessel_jn.o
	$(CXX) GWB_SMBHB.o bessel/CBESSEL.o bessel/fortran_bessel_jn.o -o $@ $(LIBS) $(RPATH)

GWB_SMBHB.o: GWB_SMBHB.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

bessel/CBESSEL.o: bessel/toms644/CBESSEL.F90
	$(FC) $(FFLAGS) -c $< -o $@

bessel/fortran_bessel_jn.o: bessel/fortran_bessel_jn.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o bessel/*.o GWB_SMBHB
