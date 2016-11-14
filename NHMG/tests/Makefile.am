########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2015
########################################################################
LIBNHMGDIR = ../src

SUFFIXES = .F90 .f90 .o .mod

.f90.mod:
	$(FC) $(FCFLAGS) -c $<

.f90.o:
	$(FC) $(FCFLAGS) -c -o $@ $<

.c.o:
	$(CC) $(CFLAGS) -c  $@ $<

clean-local:
	-rm -f *.mod *.o

cleanf:
	-rm -f *.nc fort.* *.dump *.xml gmon.out *.log *__genmod*

##
SUBDIRS	       = .
examplesdir    = $(prefix)
examples_DATA  = nh_namelist
EXTRA_DIST     = $(examples_DATA)

bin_PROGRAMS   = \
	testcuc      \
    testseamount

##
##
testcuc_SOURCES =           \
    $(LIBNHMGDIR)/libnhmg.a \
	nhydro.f90              \
	mg_testcuc.f90
##
testseamount_SOURCES=       \
	$(LIBNHMGDIR)/libnhmg.a \
	nhydro.f90              \
	mg_testseamount.f90
