# ================================
#       COMPILER & OPTIONS
# ================================
FC := mpifort

# Feature toggles
DEBUG        ?= 1    # 1 = debug, 0 = release
USE_OPENACC  ?= 0
USE_OPENMP   ?= 0
USE_NVTX     ?= 1

# Debugging toggles
ACC_DEBUG    ?= 0    # OpenACC debugging (GNU-supported only)
OMP_DEBUG    ?= 0    # OpenMP debugging

# ================================
#       BASE FLAGS
# ================================
CPPFLAGS = -I. -I$(NETCDF_FORTRAN_INC)
FCFLAGS  = -fpic -cpp
LDFLAGS  = -L$(NETCDF_FORTRAN_LIB) -lnetcdff -lnetcdf

# ================================
#       DEBUG / RELEASE FLAGS
# ================================
DEVFLAGS = -O0 -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none
RELFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops

ifeq ($(DEBUG),1)
  FCFLAGS += $(DEVFLAGS)
else
  FCFLAGS += $(RELFLAGS)
endif

# ================================
#       OPENACC SUPPORT (GNU)
# ================================
ifeq ($(USE_OPENACC),1)
  FCFLAGS  += -fopenacc
  CPPFLAGS += -D_OACC
LDFLAGS  = -L$(NETCDF_FORTRAN_LIB) -lnetcdff -lnetcdf -lm
  ifeq ($(ACC_DEBUG),1)
    FCFLAGS += -O0 -g -fcheck=all -fbacktrace \
               -fopt-info-all -fopt-info-optimized-omp
    # For runtime debug:
    # export GOMP_DEBUG=1
  endif
endif

# ================================
#       OPENMP SUPPORT
# ================================
ifeq ($(USE_OPENMP),1)
  FCFLAGS  += -fopenmp
  CPPFLAGS += -D_OMP

  ifeq ($(OMP_DEBUG),1)
    FCFLAGS += -fcheck=bounds -fcheck=pointer \
               -fsanitize=address -fsanitize=undefined \
               -fopt-info-all -fopt-info-optimized-omp
  endif
endif

# ================================
#       NVTX SUPPORT
# ================================
ifeq ($(USE_NVTX),1)
  CPPFLAGS += -DUSE_NVTX -I$(CUDA_HOME)/include
  LDFLAGS  += -L$(CUDA_HOME)/lib64 -lnvToolsExt
endif

# ================================
#       BUILD RULES
# ================================
all: model

# Compile modules in correct order
parallel_timer.o: parallel_timer.f90
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_parameters.o: module_parameters.f90
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_nvtx.o: module_nvtx.F90
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_types.o: module_types.F90 module_parameters.o
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_physics.o: module_physics.f90 parallel_timer.o module_types.o
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

module_output.o: module_output.F90 module_physics.o
        $(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

# Link final executable
model: model.f90 parallel_timer.o module_nvtx.o module_parameters.o module_types.o module_physics.o module_output.o
        $(FC) $(FCFLAGS) $(CPPFLAGS) $< -o $@ *.o $(LDFLAGS)

# Clean
clean:
        rm -f *.o *.mod model *.nc

~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
~                                                         