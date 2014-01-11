# Automake file for Simpledataio
# Written by:
# 	Edmund Highcock (edmund.highcock@users.sourceforge.net)

ACLOCAL_AMFLAGS=-I m4
# C and Fortran Libraries
lib_LTLIBRARIES = libsimpledataio.la libsimpledataiof.la

# If the compiler generates .mod files
if HAS_F90_MODULES
# If the .mod files have uppercase names
if UPPER_CASE_MOD
SDAT_MOD=SIMPLEDATAIO.mod
SDAT_WRITE_MOD=SIMPLEDATAIO_WRITE.mod
else
SDAT_MOD=simpledataio.mod
SDAT_WRITE_MOD=simpledataio_write.mod
endif
# The Cray compiler, for example, stores the module information
# in the object files, hence no modules are necessary
else
SDAT_MOD =
SDAT_WRITE_MOD =
endif
# Public headers/fortran modules
include_HEADERS = include/simpledataio.h $(SDAT_MOD) $(SDAT_WRITE_MOD)
# Sources required
libsimpledataio_la_SOURCES = src/simpledataio.c 
libsimpledataiof_la_SOURCES = src/simpledataiof.F90 src/simpledataio_write.f90

# Test programs
bin_PROGRAMS=$(top_builddir)/test/testprog 
bin_PROGRAMS+=$(top_builddir)/test/testprog_fortran

# Test program sources
__top_builddir__test_testprog_SOURCES=test/test.c
__top_builddir__test_testprog_LDADD=libsimpledataio.la
__top_builddir__test_testprog_fortran_SOURCES=test/testf.f90
__top_builddir__test_testprog_fortran_LDADD=libsimpledataiof.la libsimpledataio.la 

# Weird that this is necessary, but oh well
AM_FCFLAGS=-I/usr/include

# To get the module files
AM_FCFLAGS+=-Iinclude/

# On HPC systems with Intel/PGI compilers, force static linking
if FORCE_STATIC
AM_LDFLAGS=-all-static
else
endif

# Serial tests for C and Fortran... invoked by make check
TESTS = test/testprog test/testprog_fortran
TESTS_ENVIRONMENT = 

# If the user has specified --enable-parallel
if PARALLEL
AM_CPPFLAGS=-DPARALLEL
bin_PROGRAMS+=$(top_builddir)/test/testprog_fortran_parallel
TESTS +=  test/test_parallel.sh

__top_builddir__test_testprog_fortran_parallel_SOURCES=test/test_parallelf.f90
__top_builddir__test_testprog_fortran_parallel_LDADD=libsimpledataiof.la libsimpledataio.la
#__top_builddir__test_testprog_fortran_parallel_LDFLAGS=-lmpi -lmpi_f90
endif

CLEANFILES=test.cdf testfile.cdf test_parallel.cdf $(SDAT_MOD) $(SDAT_WRITE_MOD) include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)


#TEST_EXTENSIONS = .sh
#SH_LOG_COMPILER = bash
#LOG_COMPILER = sh

if HAS_F90_MODULES
$(SDAT_MOD): simpledataiof.lo
include/$(SDAT_MOD): $(SDAT_MOD)
	cp $(SDAT_MOD) include/.
$(SDAT_WRITE_MOD): simpledataio_write.lo
include/$(SDAT_WRITE_MOD): $(SDAT_WRITE_MOD)
	cp $(SDAT_WRITE_MOD) include/.
#simpledataio_write.o: include/$(SDAT_MOD) 
test/testf.f90: include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)
test/test_parallelf.f90: include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)
endif

# Need to generate the module file for simpledataiof before building simpledataio_write
src/simpledataio_write.f90: src/generate_simpledataio_write.rb $(SDAT_MOD)
	env ruby $< > $@

