# Location of RTE+RRTMGP libraries, module files.
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
-include $(RRTMGP_BUILD)/Makefile.conf

# Location of netcdf C and Fortran libraries. Could specify with environment variables if file doesn't exist
-include $(RRTMGP_ROOT)/examples/rfmip-clear-sky/Makefile.libs
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

#
# Build with GPTL timing if GPTL_ROOT is set
#
ifdef GPTL_ROOT
    FCFLAGS += -DUSE_TIMING -I$(GPTL_ROOT)/include
    LDFLAGS += -L$(GPTL_ROOT)/lib 
    LIBS += -lgptl -lgptl
endif

VPATH = ./:$(RRTMGP_BUILD):$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/extensions/cloud_optics:$(RRTMGP_ROOT)/examples

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%.o: %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.o
ADDITIONS += mo_garand_atmos_io.o

#
# Targets
#
all: rrtmgp_clouds

rrtmgp_clouds: $(ADDITIONS) rrtmgp_clouds.o

rrtmgp_clouds.o: $(ADDITIONS) rrtmgp_clouds.F90

mo_cloud_optics.o: mo_cloud_optics.F90
mo_load_coefficients.o: mo_simple_netcdf.o mo_load_coefficients.F90
mo_garand_atmos_io.o:   mo_simple_netcdf.o mo_garand_atmos_io.F90
mo_load_cloud_coefficients.o: mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.F90

clean:
	-rm *.o *.mod
