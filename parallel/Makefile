# ================================
#       COMPILER & OPTIONS
# ================================
# The compiler
F_COMP					:= $(MPIF90)
CPP_COMP				:= $(MPICXX) -std=c++17

# Feature toggles
DEBUG        			?= 0    # 1 = debug, 0 = release
USE_OPENACC  			?= 0
USE_OPENMP   			?= 1

# GPU architecture
CUDA_ARCH				= -gpu=cc80


# ================================
#       BASE FLAGS
# ================================
CPPFLAGS 				= -I. -I$(NETCDF_FORTRAN_INC)
FCFLAGS  				= -fpic -cpp -I$(CUDA_HOME)/include
LDFLAGS  				= -L$(NETCDF_FORTRAN_LIB) -lnetcdff -lnetcdf -lstdc++ -lm -L$(CUDA_HOME)/lib64 -lnvToolsExt


# ================================
#       DEBUG / RELEASE FLAGS
# ================================
DEVFLAGS 				= -O0 -g -C -Wall -Wextra
RELFLAGS 				= -O3 -fast -tp=native


ifeq ($(DEBUG),1)
	FCFLAGS += $(DEVFLAGS)

	ifeq ($(USE_OPENACC),1)
		FCFLAGS += -gpu=debug -traceback -Mdclchk
	endif
else
	FCFLAGS += $(RELFLAGS)
endif

# ================================
#       OPENACC SUPPORT (GNU)
# ================================
ifeq ($(USE_OPENACC),1)
	FCFLAGS  += -acc $(CUDA_ARCH) -Minfo=accel
	CPPFLAGS += -D_OACC
endif


# ================================
#       OPENMP SUPPORT
# ================================
ifeq ($(USE_OPENMP),1)
	FCFLAGS  += -MP
	CPPFLAGS += -D_OMP
endif


# ================================
#       BUILD RULES
# ================================
all: model

# Do the final linking
model: model.f90  module_parameters.o module_types.o module_physics.o module_output.o parallel_timer.o parallel_timer_c.o module_nvtx.o
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules in correct order
module_nvtx.o: module_nvtx.F90
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

parallel_timer.o: parallel_timer.f90 parallel_timer_c.o
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

parallel_timer_c.o : parallel_timer_c.cpp
	$(CPP_COMP) -c $<

module_parameters.o: module_parameters.f90
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

module_types.o: module_types.F90 module_parameters.o parallel_timer.o module_nvtx.o
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

module_physics.o: module_physics.f90 module_types.o module_nvtx.o
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

module_output.o: module_output.F90 module_physics.o module_nvtx.o
	$(F_COMP) $(FCFLAGS) $(CPPFLAGS) -c $<

# Clean
clean:
	rm -f *.o *.mod model *.nc *.txt