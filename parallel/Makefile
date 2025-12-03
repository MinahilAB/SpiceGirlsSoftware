# ================================
#       COMPILER & OPTIONS
# ================================
FC						= mpifort

# Feature toggles
DEBUG        			?= 1    # 1 = debug, 0 = release
USE_OPENACC  			?= 0
USE_OPENMP   			?= 0
USE_NVTX     			?= 1


# ================================
#       BASE FLAGS
# ================================
CPPFLAGS 				= -I. -I$(NETCDF_FORTRAN_INC)
FCFLAGS  				= -fpic -cpp
LDFLAGS  				= -L$(NETCDF_FORTRAN_LIB) -lnetcdff -lnetcdf -lstdc++ -lm


# ================================
#       DEBUG / RELEASE FLAGS
# ================================
DEVFLAGS 				= -O0 -g -Wall -Wextra -fimplicit-none -fcheck=all -fbacktrace
RELFLAGS 				= -O3 -march=native -mtune=native -Ofast -ffast-math -funroll-loops

ifeq ($(DEBUG),1)
	FCFLAGS += $(DEVFLAGS)
	
	ifeq ($(USE_OPENACC),1)
		FCFLAGS += -fopt-info-all
		FCFLAGS  := $(filter-out -fcheck=all -fbacktrace,$(FCFLAGS))
	else ifeq ($(USE_OPENMP),1)
		FCFLAGS += -fopt-info-optimized-omp
	endif
else
	FCFLAGS += $(RELFLAGS)
endif


# ================================
#       OPENACC SUPPORT (GNU)
# ================================
ifeq ($(USE_OPENACC),1)
	FCFLAGS  += -fopenacc -foffload=nvptx-none=-lm
	CPPFLAGS += -D_OACC
endif


# ================================
#       OPENMP SUPPORT
# ================================
ifeq ($(USE_OPENMP),1)
	FCFLAGS  += -fopenmp
	CPPFLAGS += -D_OMP
endif


# ================================
#       NVTX SUPPORT
# ================================
ifeq ($(USE_NVTX),1)
	FCFLAGS += -I$(CUDA_HOME)/include
	CPPFLAGS += -DUSE_NVTX
	LDFLAGS  += -L$(CUDA_HOME)/lib64 -lnvToolsExt
endif


# ================================
#       BUILD RULES
# ================================

SRC      := $(wildcard *.f90 *.F90)
OBJ      := $(SRC:.f90=.o)
OBJ      := $(OBJ:.F90=.o)

all: model

# Do the final linking
model: $(OBJ)
	$(FC) $(FCFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules in correct order

ifeq ($(USE_NVTX),1)
module_nvtx.o: module_nvtx.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<
endif

parallel_timer.o: parallel_timer.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_parameters.o: module_parameters.f90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_types.o: module_types.F90 module_parameters.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_physics.o: module_physics.f90 module_types.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_output.o: module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

# Clean
clean:
	rm -f *.o *.mod model *.nc *.txt