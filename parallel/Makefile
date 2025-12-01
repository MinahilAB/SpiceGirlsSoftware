FC					:= mpifort

# Build configuration toggles
DEBUG       		?= 1
USE_OPENACC 		?= 0
USE_OPENMP  		?= 0
USE_NVTX            ?= 1

# Development/release flags
# 
DEVFLAGS 			:= -O0 -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none
RELFLAGS			:= -O3 -march=native -mtune=native -ffast -ffast-math -funroll-loops

# Base flags
CPPFLAGS			= -I. -I$(NETCDF_FORTRAN_INC)
FCFLAGS				= -fpic -cpp
LDFLAGS				= -L$(NETCDF_FORTRAN_LIB) -lnetcdff -lnetcdf

# Apply flags according to specified options
ifeq ($(DEBUG),1)
  FCFLAGS += $(DEVFLAGS)
else
  FCFLAGS += $(RELFLAGS)
endif

# OpenACC support
ifeq ($(USE_OPENACC),1)
  FCFLAGS  += -fopenacc
  CPPFLAGS += -D_OACC
endif

# OpenMP support
ifeq ($(USE_OPENMP),1)
  FCFLAGS  += -fopenmp
  CPPFLAGS += -D_OMP
endif

# NVTX support
ifeq ($(USE_NVTX),1)
  CPPFLAGS += -DUSE_NVTX -I$(CUDA_HOME)/include
  LDFLAGS  += -L$(CUDA_HOME)/lib64 -lnvToolsExt
endif

# Build rules
all: model

model : model.f90 module_nvtx.F90 module_nvtx.o module_parameters.o module_types.F90 module_physics.o module_output.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -o $@ *.o $(LDFLAGS)

module_nvtx.o : module_nvtx.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -c

module_parameters.o : module_parameters.f90
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -c

module_types.o : module_types.F90 module_parameters.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -c

module_output.o : module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -c

module_physics.o : module_physics.f90 module_types.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) $< -c

clean :
	rm -f *.o *.mod model *.nc
