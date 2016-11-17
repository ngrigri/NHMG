########################################################################
##- Author@mail: Nicolas.Grima@univ-brest.fr
##- Date       : November 2016
########################################################################
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

SRCDIR = ../src
LIBNHMGDIR = ${SRCDIR}

##
SUBDIRS	       = .
testsdir    = $(prefix)
tests_DATA  = nhmg_namelist
EXTRA_DIST     = $(tests_DATA)

bin_PROGRAMS   = \
	testcuc      \
    testseamount

##
##
testcuc_SOURCES =  \
	mg_testcuc.f90 

testcuc_LDADD = $(LIBNHMGDIR)/libnhmg.a

##
testseamount_SOURCES=   \
	mg_testseamount.f90

testseamount_LDADD = $(LIBNHMGDIR)/libnhmg.a
