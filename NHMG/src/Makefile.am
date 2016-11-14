########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2015
########################################################################
FORTRAN_TOOLS=.

SUFFIXES = .F90 .f90 .o .mod

.f90.mod:
	$(FC) -c $<

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
	mg_mpi.f90             \
	mg_tictoc.f90          \
	mg_namelist.f90        \
	mg_netcdf_out.f90      \
	mg_grids.f90           \
	mg_mpi_exchange.f90    \
	mg_gather.f90          \
    mg_zr_zw.f90           \
	mg_define_matrices.f90 \
	mg_compute_rhs.f90     \
    mg_correct_uvw.f90     \
    mg_set_bbc.f90         \
    mg_compute_barofrc.f90 \
    mg_compute_fluxes.f90  \
	mg_relax.f90           \
	mg_intergrids.f90      \
    mg_btbc_coupling.f90   \
	mg_solvers.f90
