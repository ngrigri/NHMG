########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2016
########################################################################
SRCDIR = ../src
LIBNHMGDIR = ${SRCDIR}

AM_FCFLAGS  = ${AM_CPPFLAGS} $(FCFLAGS_f90) -I$(top_srcdir)/src
AM_FFLAGS   = ${AM_FCFLAGS}
##AM_LDFLAGS  = -L/opt/linux/mpich2/1.5-intel12.1/lib -lmpich 
##AM_LDFLAGS += -L/opt/linux/netcdf/4.1.3-intel12/64/lib -lnetcdf -lnetcdff
AM_LDFLAGS  = $(TESTLIBS)
AM_LDFLAGS += ${top_builddir}/src/libnhmg.la -lnhmg

F90_TESTS = testseamount testcuc

check_PROGRAMS = $(F90_TESTS)
TESTS          = $(F90_TESTS)


##
SUBDIRS	       = .
testsdir    = $(prefix)
tests_DATA  = nhmg_namelist ts_namelist cuc_namelist
EXTRA_DIST     = $(tests_DATA)

bin_PROGRAMS   = \
    testseamount \
	testcuc

##
testseamount_SOURCES=       \
	mg_zr_zw.f90            \
	mg_mpi_exchange_ijk.f90 \
    mg_setup_tests.f90      \
	mg_testseamount.f90

##testseamount_LDADD = $(LIBNHMGDIR)/libnhmg.a
##
testcuc_SOURCES =           \
	mg_zr_zw.f90            \
	mg_mpi_exchange_ijk.f90 \
	mg_setup_tests.f90      \
	mg_testcuc.f90 

##testcuc_LDADD = $(LIBNHMGDIR)/libnhmg.a

clean-local:
	-rm -f *.mod *.o

cleanf:
	-rm -f *.nc fort.* *.dump *.xml gmon.out *.log *__genmod*

##
