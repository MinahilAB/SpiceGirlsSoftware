FC = gfortran
#FCFLAGS = -O0 -Wall -pedantic -fcheck=all -fbacktrace -I. `nf-config --fflags`
FCFLAGS = -Ofast -I. `nf-config --fflags` -fPIC
LDFLAGS = -L`nc-config --libdir` `nf-config --flibs`

all: model

model : model.f90 module_parameters.o module_types.F90 module_physics.o \
        module_output.o
	$(FC) $(FCFLAGS) $< -o $@ *.o $(LDFLAGS)

module_parameters.o : module_parameters.f90
	$(FC) $(FCFLAGS) $< -c

module_types.o : module_types.F90 module_parameters.o
	$(FC) $(FCFLAGS) $< -c

module_output.o : module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) $< -c

module_physics.o : module_physics.f90 module_types.o
	$(FC) $(FCFLAGS) $< -c

clean :
	rm -f *.o *.mod model *.nc
