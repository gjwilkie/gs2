RANLIB=ranlib
# If CPU is undefined; use the test_os script to figure it out:
ifeq ($(CPU),)
   CPU_tmp := $(shell sh ./test_os)
   ifeq ($(CPU_tmp),ALPHA)
      CPU := ALPHA
#     CPU := ALPHA_NAG
   else
      CPU := $(CPU_tmp)
   endif
endif

DATE=$(shell date +%D | sed 's/\///g')

##################################################### PLATFORM DEPENDENT FLAGS

# (in gnu make, $$ translates to $, and shell then does variable substitution)

# T3E options:
ifeq ($(CPU),T3E)
  FC = f90
  PLATFORM_LINKS = t3e
  FFLAGS =  -N80 -I/usr/local/include -M1110 
  F90FLAGS =  -I/usr/local/include -M1110 
  FLIBS = -lmpi -L/usr/local/lib  $$NETCDF 

  ifneq ($(debug),on)
    FFLAGS += -O vector3 
    F90FLAGS +=  -O vector3
  else
    FFLAGS   += -g -R abcs -e i
    F90FLAGS += -g -R abcs -e i
    FLIBS   += -Wl"-D preset=inf" # -lmalloc 
  endif

endif

# Cray Fortran Compiler Options:
# -Rabcs, run time checks for arguments, bounds, array conformance, strings
# -N80 accepts 80 character lines
# -g debug with no optimization
# -eA (and link with -lapp) to use apprentice
# -ei check for some kinds of uninitialized variables
# -e0 initialize local stack vars to zero
# -M1110 ignore warning messages about treating "double precision" as
#        "only" 64 bits
#
# when running with code linked with -lmalloc, do:
# setenv MEMCHK 1  # check heap correctness after every memory utility call
# setenv MEMINDEF 1  # initialize malloc memory to NAN and invalid pointers.

# C90/J90 options:
ifeq ($(CPU),C90)
  FC = f90
  PLATFORM_LINKS = c90
  FFLAGS =  -N80 -I/usr/local/include -M1110 
  F90FLAGS =  -I/usr/local/include -M1110
  FLIBS = -L/usr/local/lib $$NETCDF

  ifneq ($(debug),on)
    FFLAGS +=  -O vector3 
    F90FLAGS +=  -O vector3  
  else
    FFLAGS   += -g -R abcs -e i
    F90FLAGS += -g -R abcs -e i
    FLIBS   += -Wl"-D preset=inf" # -lmalloc
  endif

endif

# J90 options so that variables are private to each process (needed for MPI):
# (I can't get this to work.  So for now drop the "-a taskcommon" switch,
# which restricts us to 1 processor on the J90.)
#  FFLAGS =  -N80 -I/usr/local/include -M1110 -O vector3 -a taskcommon
#  F90FLAGS =  -I/usr/local/include -M1110 -O vector3 -a taskcommon
#  FLIBS_itg = -L/usr/local/lib -lnag $$NETCDF
#
# debug options:
#  FFLAGS =  -N80 -I/usr/local/include -M1110 -g -R abcs -e i
#  F90FLAGS =     -I/usr/local/include -M1110 -g -R abcs -e i
#  FLIBS_itg = -L/usr/local/lib -lnag $$NETCDF -lmalloc -Wl"-D preset=inf"
#                 # -lapp

# SGI Origin-2000 options:
ifeq ($(CPU),SGI)
  FC = f90
  FLIBS  =    -L/usr/pppl/lib -lnetcdf -lmpi -lscs
  PLATFORM_LINKS = origin
  FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 
  F90FLAGS =  -I/usr/pppl/include -64 -r8 

  ifneq ($(debug),on)
    FFLAGS +=   -O -TARG:platform=ip27
    F90FLAGS += -O -TARG:platform=ip27
# Other options tried (no more than 15% speedup):
#    FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
#    F90FLAGS =  -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
  else
    FFLAGS +=  -g -DEBUG:div_check=3:trap_uninitialized=on
    F90FLAGS +=  -g -DEBUG:div_check=3:trap_uninitialized=on
#       -DEBUG:div_check=3:subscript_check=on:trap_uninitialized=on
# should be the full debug options, but then the compiler bombs on some
# of the routines (even the new beta version 7.3 of the compiler).
  endif

endif

# NERSC IBM options:
ifeq ($(CPU),RS6000)
  FC = mpxlf90_r
  PLATFORM_LINKS = ibm
  FFLAGS   = -qrealsize=8
  F90FLAGS = -qrealsize=8 -qsuffix=f=f90
  FLIBS = $$NETCDF 
  ifneq ($(debug),on)
#    FFLAGS   += -O4
#    F90FLAGS += -O4
    FFLAGS   += -O3 -qarch=pwr3 -qtune=pwr3
    F90FLAGS += -O3 -qarch=pwr3 -qtune=pwr3
  else
    FFLAGS   += -g 
    F90FLAGS += -g 
    FLIBS    += # $$TRACE_MPIF
  endif
endif

# DEC alpha options:
ifeq ($(CPU),ALPHA)
  FC = f95
  FLIBS = -L/usr/local/lib -lnetcdf -L/usr/lib
# FLIBS  = -L/usr/local/lib -lnagdx -lnetcdf -ldxml
  PLATFORM_LINKS = alpha
  FFLAGS   = -r8 -extend_source
  F90FLAGS = -r8 

  ifeq ($(debug),on)
     FFLAGS   += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -check output_conversion
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -check output_conversion
  else
     FFLAGS   +=  -O -fast -w 
     F90FLAGS +=  -O -fast -w
  endif

endif

# options for Linux on a DEC alpha with DEC/Compaq F90:
ifeq ($(CPU),LINUX_alpha)
  FC = f90
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf
  PLATFORM_LINKS = linux_alpha
  FFLAGS   = -r8 -extend_source
  F90FLAGS = -r8 

  ifeq ($(debug),on)
     FFLAGS   += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -align sequence
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -align sequence
  else
     FFLAGS   +=  -O -fast -w 
     F90FLAGS +=  -O -fast -w
  endif

endif

# options for NAG f95 on a DEC alpha
ifeq ($(CPU),ALPHA_NAG)
  FC = /usr/local/bin/f95
  FLIBS = -L/usr/local/lib -lnetcdf
  PLATFORM_LINKS = alpha_nag
  FFLAGS = -132 -I/usr/local/include -r8 
  F90FLAGS =    -I/usr/local/include -r8 
  ifeq ($(debug),on)
    FFLAGS   += -C -g90 -dusty
    F90FLAGS += -C -g90 -dusty
  else
    FFLAGS   += -O
    F90FLAGS += -O 
  endif
endif

# options for Linux with pgf90:
ifeq ($(CPU),LINUX_pg)
  FC = pgf90
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf
  PLATFORM_LINKS = linux
  FFLAGS =  -r8 -132 -I/usr/local/include -Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp
  F90FLAGS = -r8 -I/usr/local/include -module ../mod -Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp

  ifeq ($(debug),on)
    FFLAGS   += -g -Mbounds
    F90FLAGS += -g -Mbounds
  else
    FFLAGS   += -w -O 
    F90FLAGS += -w -O
  endif

endif

# options for Linux with Fujitsu f90
ifeq ($(CPU),LINUX_fuj)
  FC = mpif90
  FLIBS = -L/usr/local/lib -lnetcdf
  PLATFORM_LINKS = linux_fuj
  FFLAGS   = -w   -C cdRR8 -I/usr/local/include -X9 -static-flib -Kfast 
  F90FLAGS = -A m -C cdRR8 -I/usr/local/include -X9 -static-flib -Kfast 
  ifeq ($(debug),on) 
    FFLAGS   += -g -H easu
    F90FLAGS += -g -H easu
  else
    FFLAGS   += -O -f 2004,2006,2008 -Wa,--no-warn
    F90FLAGS += -O -f 2004,2006,2008 -Wa,--no-warn
  endif
endif

# options for Linux with NAG f95:
ifeq ($(CPU),LINUX)
  FC = f95
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf
  PLATFORM_LINKS = linux
  FFLAGS   = -r8 -132 -w -f77 -I/usr/local/include
  F90FLAGS = -r8 -w -f77 -I/usr/local/include -mismatch

  ifeq ($(debug),on)
    FFLAGS   += -C -g90 -gline 
    F90FLAGS += -C -g90 -gline 
  else
    FFLAGS   += -O 
    F90FLAGS += -O
  endif

endif

util_obj = spl.o mds.o prec.o netcdf.o 

ifeq ($(CPU),T3E)
  util_obj += mpptime.o
endif

# always build the dummy libraries mdslib.a and netcdf_stub.a
# to have available as options:
all:	utils.a mdslib.a 

utils.a: $(util_obj)
	$(AR) rc $@ $(util_obj)
	$(RANLIB) $@

mds.o: prec.o

mdslib.a: mdslib.o
	$(AR) rc $@ mdslib.o
	$(RANLIB) $@

clean:
	rm -f *.o *.mod *.a

mpptime.o: 
	$(CC) -c mpptime.c

######################################################################## RULES
.SUFFIXES:
.SUFFIXES: .f90 .f

.f90.o: 
	$(FC) $(F90FLAGS) -c $<

.f.o: 
	$(FC) $(FFLAGS) -c $<

################################################################### DIRECTIVES
# search path where gnu make will look for files:
# VPATH = ./src


# If no other rules are found, use the defaults:

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

test_make:
	@echo FC is $(FC)
	@echo FFLAGS is $(FFLAGS)
	@echo F90FLAGS is $(F90FLAGS)
	@echo debug is $(debug)
	@echo CPU is $(CPU)

