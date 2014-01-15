
# Automake file for Simpledataio
# Written by:
# 	Edmund Highcock (edmund.highcock@users.sourceforge.net)

DEFS=

#$(info "DEFS is (DEFS)")
MYDEFS=@MYDEFS@
MYDEFSF=@MYDEFSF@
#CDEFS=$(MYDEFS)
#$(info "DEFS is $(DEFS)")
#DEFS=
AM_CFLAGS=$(MYDEFS)
#FDEFS=$(subst -D,$(FDEF),$(MYDEFS))
AM_FCFLAGS=$(MYDEFSF)
ACLOCAL_AMFLAGS=-I m4

# Flag needed for defining macros to the fortran preprocessor (usually -D)
FDEF=@FC_DEFINE@

# Now adjust fortran and C preprocessor flags separately
#
#
if !ARCHIVEONLY
LA = la
LO = lo
else
LA = a
LO = o
endif

# C and Fortran Libraries
if !ARCHIVEONLY
lib_LTLIBRARIES = libsimpledataio.la libsimpledataiof.la
lib_LIBRARIES =
else
lib_LIBRARIES = libsimpledataio.a libsimpledataiof.a
#lib_LTLIBRARIES=
endif

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
if !ARCHIVEONLY
libsimpledataio_la_SOURCES = src/simpledataio.c
libsimpledataiof_la_SOURCES = src/simpledataiof.F90 src/simpledataio_write.F90
else
libsimpledataio_a_SOURCES = src/simpledataio2.c 
libsimpledataiof_a_SOURCES = src/simpledataiof2.F90 src/simpledataio_write2.F90
endif

if !ARCHIVEONLY
# Test programs
bin_PROGRAMS=$(top_builddir)/test/testprog 
bin_PROGRAMS+=$(top_builddir)/test/testprog_fortran

# Test program sources
__top_builddir__test_testprog_SOURCES=test/test.c
__top_builddir__test_testprog_LDADD=libsimpledataio.$(LA)
__top_builddir__test_testprog_fortran_SOURCES=test/testf.f90
__top_builddir__test_testprog_fortran_LDADD=libsimpledataiof.$(LA) libsimpledataio.$(LA) 
endif


# To get the module files
AM_FCFLAGS+=-Iinclude/


# On HPC systems with Intel/PGI compilers, force static linking
if FORCE_STATIC
AM_LDFLAGS=-all-static
else
endif

if HAS_FORTRAN_NETCDF
AM_FCFLAGS+= $(FDEF)FORTRAN_NETCDF
endif
if HAS_ISO_C_BINDING
AM_FCFLAGS+= $(FDEF)ISO_C_BINDING
endif

# Serial tests for C and Fortran... invoked by make check
TESTS = test/testprog test/testprog_fortran
TESTS_ENVIRONMENT = 

# If the user has specified --enable-parallel
if PARALLEL
AM_CPPFLAGS=-DPARALLEL
if !ARCHIVEONLY
bin_PROGRAMS+=$(top_builddir)/test/testprog_fortran_parallel
TESTS +=  test/test_parallel.sh

__top_builddir__test_testprog_fortran_parallel_SOURCES=test/test_parallelf.f90
__top_builddir__test_testprog_fortran_parallel_LDADD=libsimpledataiof.$(LA) libsimpledataio.$(LA)
#__top_builddir__test_testprog_fortran_parallel_LDFLAGS=-lmpi -lmpi_f90
endif
endif

CLEANFILES=test.cdf testfile.cdf test_parallel.cdf $(SDAT_MOD) $(SDAT_WRITE_MOD) include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)
CLEANFILES+= src/simpledataio2.c src/simpledataiof2.F90 src/simpledataio_write2.F90



#TEST_EXTENSIONS = .sh
#SH_LOG_COMPILER = bash
#LOG_COMPILER = sh

if ARCHIVEONLY
  TWO = 2
endif
if HAS_F90_MODULES
$(SDAT_MOD): simpledataiof$(TWO).$(LO)
include/$(SDAT_MOD): $(SDAT_MOD)
	cp $(SDAT_MOD) include/.
$(SDAT_WRITE_MOD): simpledataio_write$(TWO).$(LO)
include/$(SDAT_WRITE_MOD): $(SDAT_WRITE_MOD)
	cp $(SDAT_WRITE_MOD) include/.
#simpledataio_write.o: include/$(SDAT_MOD) 
test/testf.f90: include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)
test/test_parallelf.f90: include/$(SDAT_MOD) include/$(SDAT_WRITE_MOD)
endif

if ARCHIVEONLY
src/simpledataio2.c: src/simpledataio.c
	cp $< $@
src/simpledataiof2.F90: src/simpledataiof.F90
	cp $< $@
src/simpledataio_write2.F90: src/simpledataio_write.F90
	cp $< $@
endif



# Need to generate the module file for simpledataiof before building simpledataio_write
src/simpledataio_write.F90: src/generate_simpledataio_write.rb $(SDAT_MOD)
	env ruby $< > $@


