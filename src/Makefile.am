########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2015
########################################################################
FORTRAN_TOOLS=.

SUFFIXES = .F90 .f90 .o .mod

.f90.mod:
	$(FC) $(FCFLAGS) -c $<

.f90.o:
	$(FC) -c -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c  $@ $<

clean-local:
	-rm -f *.mod *.o *.a

cleanf:
	-rm -f *.nc fort.* *.dump *.xml gmon.out *.log *__genmod*

##
SUBDIRS	       = .

lib_LTLIBRARIES = libnhmg.la
libnhmg_la_LDFLAGS = -version-number 0:0:0
libnhmg_la_SOURCES =       \
	mg_cst.f90             \
	mg_mpi.f90             \
	mg_tictoc.f90          \
	mg_namelist.f90        \
	mg_netcdf_out.f90      \
	mg_grids.f90           \
	mg_mpi_exchange.f90    \
	mg_gather.f90          \
	mg_horiz_grids.f90     \
	mg_vert_grids.f90     \
	mg_bbc.f90  \
	mg_fluxes.f90  \
	mg_projection.f90 \
	mg_relax.f90           \
	mg_intergrids.f90      \
	mg_coupling.f90   \
	mg_solvers.f90         \
	nhmg.f90

#include_HEADERS = nhmg.mod
nodist_include_HEADERS =   \
	mg_cst.mod             \
	mg_mpi.mod             \
	mg_tictoc.mod          \
	mg_namelist.mod        \
	mg_netcdf_out.mod      \
	mg_grids.mod           \
	mg_mpi_exchange.mod    \
	mg_gather.mod          \
	mg_horiz_grids.mod     \
	mg_vert_grids.mod     \
	mg_bbc.mod  \
	mg_fluxes.mod  \
	mg_projection.mod \
	mg_relax.mod           \
	mg_intergrids.mod      \
	mg_coupling.mod   \
	mg_solvers.mod         \
	nhmg.mod

#libnhmg_la_DEPENDENCIES = nhmg.mod mg_mpi.mod mg_grids.mod


