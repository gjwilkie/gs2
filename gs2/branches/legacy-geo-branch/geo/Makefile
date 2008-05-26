RANLIB=ranlib
NTDIR=..
UTILDIR=../utils
# If CPU is undefined; use the test_os script to figure it out:

ifeq ($(CPU),)
   CPU_tmp := $(shell sh ../test_os)
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
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o
  FFLAGS =  -N80 -I/usr/local/include -M1110,7212 
  F90FLAGS =     -I/usr/local/include -M1110,7212 -p$(UTILDIR)
  FLIBS = -lmpi -L/usr/local/lib  $$NETCDF 

  ifneq ($(debug),on)
    FFLAGS += -O vector3
    F90FLAGS += -O vector3
  else
    FFLAGS   +=  -g -R abcs -e i
    F90FLAGS +=  -g -R abcs -e i 
    FLIBS   +=  -Wl"-D preset=inf" # -lmalloc 
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

# C90/Joptions:
ifeq ($(CPU),C90)
  FC = f90
  PLATFORM_LINKS = c90
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o 
  FFLAGS =  -N80 -I/usr/local/include -M1110
  F90FLAGS =  -I/usr/local/include -M1110 
  FLIBS = -L/usr/local/lib  $$NETCDF 

  ifneq ($(debug),on)
    FFLAGS   += -O vector3 
    F90FLAGS += -O vector3  
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
  FLIBS = -L/usr/pppl/lib -lnetcdf -lmpi -lscs
  PLATFORM_LINKS = origin
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o 
  FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 
  F90FLAGS =  -I/usr/pppl/include -64 -r8 

  ifneq ($(debug),on)
    FFLAGS   += -O -TARG:platform=ip27
    F90FLAGS += -O -TARG:platform=ip27
# Other options tried (no more than 15% speedup):
#    FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
#    F90FLAGS =  -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
  else
    FFLAGS   += -g -DEBUG:div_check=3:trap_uninitialized=on
    F90FLAGS += -g -DEBUG:div_check=3:trap_uninitialized=on
#       -DEBUG:div_check=3:subscript_check=on:trap_uninitialized=on
# should be the full debug options, but then the compiler bombs on some
# of the routines (even the new beta version 7.3 of the compiler).
  endif

endif

# NERSC IBM options:
ifeq ($(CPU),RS6000)
  FC = mpxlf90_r
#  FC = h5pfc
  PLATFORM_LINKS = ibm
#  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o $(UTILDIR)/netcdf.o
  u_obj = $(UTILDIR)/utils.a
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o
  FFLAGS   = -qautodbl=dbl4 -I $(UTILDIR)
  F90FLAGS = -qautodbl=dbl4 -qsuffix=f=f90 -I $(UTILDIR)
  FLIBS = $$NETCDF 
  ifneq ($(debug),on)
#    FFLAGS   += -O4
#    F90FLAGS += -O4
    FFLAGS   += -O3 -qarch=auto -qtune=auto
    F90FLAGS += -O3 -qarch=auto -qtune=auto
  else
    FFLAGS   += -C -g 
    F90FLAGS += -C -g 
    FLIBS    += # $$TRACE_MPIF
  endif
endif

# imported from AGK by TT FEB 01 2008
# NERSC Cray XT4 (Franklin) options: GGH 07 SEP 17
ifeq ($(CPU),XT4)
  FC = ftn
  FC5 = ftn
  PLATFORM_LINKS = cxt4
  u_obj = $(UTILDIR)/utils.a
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o
  F90FLAGS = -r8 -target=linux -I$(UTILDIR)
#  FFLAGS = -r8 -target=linux -I$$NETCDF_DIR/include $$FFTW_INCLUDE_OPTS -I$(UTILDIR)
  FLIBS = -L$$NETCDF_DIR/lib -lnetcdf
  ifneq ($(debug),on)
#    F90FLAGS += -O2# -fastsse
    F90FLAGS += -fastsse
  else
    F90FLAGS += -g
  endif
endif
# import end

# Dawson G5 cluster:
ifeq ($(CPU),Dawson)
  FC = xlf95
  PLATFORM_LINKS = ibm
  u_obj = $(UTILDIR)/utils.a
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o
  FFLAGS   = -qautodbl=dbl -I $(UTILDIR)
  F90FLAGS = -qmoddir=/tmp/bdorland -I/tmp/bdorland -qautodbl=dbl4 -qsuffix=f=f90 \
	-I $(UTILDIR) 
  FLIBS = -L/u/local/apps/netcdf/lib -lnetcdf 
  ifneq ($(debug),on)
    F90FLAGS += -O3 -qarch=g5 -qtune=g5
  else
    F90FLAGS += -g 
  endif
endif

ifeq ($(CPU),LINUX_ifort_mrf)
  FC = ifort
  PLATFORM_LINKS = mrf
  FLIBS = -L/usr/local/lib -lnetcdf #-L/opt/intel/mkl/9.1/lib/32 -lvml
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o $(UTILDIR)/mdslib.o
  stubs = radstub.o vdimstub.o
  ifeq ($(debug),on)
    FFLAGS = -r8 -g -I/usr/local/include -check bounds -lowercase
    F90FLAGS = -r8 -g -I/usr/local/include/netcdf -module $(UTILDIR) -check bounds -lowercase
  else	
    FFLAGS   = -r8 -O3 -I $(UTILDIR) -lowercase -I/usr/local/include
    F90FLAGS = -r8 -O3 -I $(UTILDIR) -I/usr/local/include -lowercase
  endif
endif
 

# DEC alpha options:
ifeq ($(CPU),ALPHA)
  FC = f95
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf ../utils/mdslib.a
#	-L/usr/local/mdsplus_new/lib -lMdsLib -lMdsShr
# FLIBS  = -L/usr/local/lib -lnagdx -lnetcdf -ldxml
  PLATFORM_LINKS = alpha
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o
  FFLAGS   =  -I/usr/local/include -r8 -extend_source -I$(UTILDIR)
  F90FLAGS =  -I/usr/local/include -r8 -I$(UTILDIR)

  ifeq ($(debug),on)
     FFLAGS   += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -check output_conversion
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -check output_conversion
  else
     FFLAGS   += -O -fast -w
     F90FLAGS += -O -fast -w 
  endif

endif

# options for LINUX_alpha on a DEC alpha with DEC/Compaq F90:
# f90 is a link to "fort", for man pages do "man fort"
ifeq ($(CPU),LINUX_alpha)
  FC = f90
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf
  PLATFORM_LINKS = linux_alpha
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o 
  FFLAGS   =  -I/usr/local/include -r8 -extend_source -I$(UTILDIR)
  F90FLAGS =  -I/usr/local/include -r8 -I$(UTILDIR)

  ifeq ($(debug),on)
     FFLAGS   += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -align sequence
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -align sequence
  else
     FFLAGS   += -O -fast -w 
     F90FLAGS += -O -fast -w 
  endif

endif

# options for Linux with Lahey lf95
ifeq ($(CPU),LINUX_lf95)
  FC = mpif90
  FLIBS = -L/usr/local/lib -lnetcdf $(UTILDIR)/mdslib.a
  PLATFORM_LINKS = linux_lf95
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o 
  F90FLAGS = --dbl --ml cdecl -I$(UTILDIR)

  ifeq ($(static),on)
    F90FLAGS += --staticlink
  endif

  ifeq ($(debug),on) 
    F90FLAGS += -g  --chk aesu
  else
    F90FLAGS += -O 
  endif

endif
# options for NAG f95 on a DEC alpha
ifeq ($(CPU),ALPHA_NAG)
  FC = /usr/local/bin/f95
  FLIBS = -L/usr/local/lib -lnetcdf
  PLATFORM_LINKS = alpha_nag
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o 
  FFLAGS =  -C -132 -I/usr/local/include -r8 -dusty
  F90FLAGS =  -C -I/usr/local/include -r8 -dusty

  ifeq ($(debug),on)
    FFLAGS   += -g90 
    F90FLAGS += -g90
  else
    FFLAGS   += -O
    F90FLAGS += -O
  endif

endif

# options for Linux with pgf90:  NOT RECOMMENDED NOT RECOMMENDED
ifeq ($(CPU),LINUX_pg)
  FC = pgf90
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf 
  PLATFORM_LINKS = linux
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o $(UTILDIR)/mdslib.o
  stubs = radstub.o vdimstub.o 
  ifeq ($(debug),on)
    FFLAGS = -r8 -I/usr/local/include -g -Mbounds -Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp
    F90FLAGS = -g -r8 -I/usr/local/include -module ../mod -Mbounds \
	-Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp
  else
    FFLAGS = -w -r8 -I/usr/local/include -O \
	-Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp
    F90FLAGS = -w -r8 -module ../mod -I/usr/local/include \
	-Mnoupcase -Mdalign -Mdefaultunit -Ktrap=fp
  endif
endif

# options for Linux with Absoft f90
ifeq ($(CPU),LINUX_abs)
  FC = mpif90
  FLIBS = -L/usr/lib -lnetcdf
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o $(UTILDIR)/mdslib.o
  PLATFORM_LINKS = linux_fuj
  F90FLAGS_base = -N113 -p$(UTILDIR)

  F90FLAGS_0 = $(F90FLAGS_base) -YEXT_SFX=
  F90FLAGS_1 = $(F90FLAGS_base) -YEXT_SFX=_
  F90FLAGS_2 = $(F90FLAGS_base) -YEXT_SFX=__

  F90FLAGS = $(F90FLAGS_base)

  ifeq ($(debug),on) 
    F90FLAGS += -g -Rbcs
    F90FLAGS_0 += -g -Rbcs
    F90FLAGS_1 += -g -Rbcs
    F90FLAGS_2 += -g -Rbcs
  else
    F90FLAGS += -O 
    F90FLAGS_0 += -O 
    F90FLAGS_1 += -O 
    F90FLAGS_2 += -O 
  endif

endif

# options for Linux with Fujitsu f90
ifeq ($(CPU),LINUX_fuj)
  FC = mpif90
  FLIBS = /usr/local/lib -lnetcdf \
	../utils/mdslib.a -L/u/bdorland/fftw/lib -lfftw -lrfftw
  PLATFORM_LINKS = linux_fuj
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o
  FFLAGS   = -w   -C cdRR8 -I/usr/local/include -X9 -I$(UTILDIR) -static-flib -Kfast 
  F90FLAGS = -A m -C cdRR8 -I/usr/local/include -X9 -I$(UTILDIR) -static-flib -Kfast 

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
  FC = mpif90
  FLIBS = -L/usr/local/lib -lnetcdf $(UTILDIR)/mdslib.a
#	-L/usr/local/mdsplus/lib -lMdsLib -lMdsShr 
  PLATFORM_LINKS = linux
  u_obj = $(UTILDIR)/spl.o $(UTILDIR)/mds.o 
  stubs = radstub.o vdimstub.o 
  FFLAGS = -r8 -132 -mismatch -I/usr/local/include  
  F90FLAGS = -r8 -mismatch -I$(UTILDIR) -I/usr/local/include

  ifeq ($(debug),on)
    FFLAGS   += -C=array -C=bits -C=dangling -C=do -C=present -C=pointer -gline 
    F90FLAGS += -C=array -C=bits -C=dangling -C=do -C=present -C=pointer -gline 
  else
    FFLAGS   += -O4 -w
    F90FLAGS += -O4 -w
  endif

endif

# TT> local intel linux machine options:
ifeq ($(CPU),LINUX_ifort)
  FC = ifort
  F90FLAGS = -r8 -I$(UTILDIR)
  FLIBS = -lnetcdf -lsvml
  PLATFORM_LINKS = linux
  u_obj = $(UTILDIR)/utils.a $(UTILDIR)/mdslib.a
  stubs = radstub.o vdimstub.o 
  ifeq ($(debug),on)
    F90FLAGS += -g -implicitnone -warn all -check bounds -traceback
  else
# -xN is for Pentium 4 chips
    F90FLAGS += -O3 -xN -vec_report0
  endif
endif
# <TT

geo_obj = geometry.o geq.o veq.o eeq.o peq.o leq.o deq.o ideq.o

all:	mods geo.a eiktest ball

geo.a: $(stubs) $(geo_obj)
	$(AR) rc $@ $(geo_obj) $(stubs)
	$(RANLIB) $@

clean:
	rm -f *.o *~ geometry.a *.tmp *.mod *.trace core eiktest
distclean:
	-rm ball eiktest

configure: configure.in
	autoconf --local=config

geometry.o : veq.o geq.o eeq.o leq.o peq.o deq.o radstub.o ideq.o

ball: ball.o geo.a $(u_obj) $(stubs)
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(FFLAGS) -o ball ball.o $(stubs) $(FLIBS) ;; \
		*)  $(FC) $(FFLAGS) -o ball ball.o $(stubs) $(geo_obj) $(u_obj) $(FLIBS) ;; \
	esac

eiktest: et.o geo.a $(u_obj) $(stubs)
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) -o eiktest et.o $(stubs) $(FLIBS) ;; \
		*) $(FC) $(FFLAGS) -o eiktest et.o $(stubs) $(geo_obj) $(u_obj) $(FLIBS) ;; \
	esac 

ts : testsmoothf.o smoothf.o
	$(FC) $(FFLAGS) -o ts testsmoothf.o smoothf.o 

ifeq ($(CPU),LINUX_abs)

peq.o: 
	$(FC) $(F90FLAGS_1) -c peq.f90
geq.o: 
	$(FC) $(F90FLAGS_1) -c geq.f90

endif

et.o: geometry.o 
ball.o: geometry.o 

peq.o geq.o: $(UTILDIR)/utils.a
veq.o: $(UTILDIR)/utils.a vdimstub.o
eeq.o: $(UTILDIR)/utils.a
deq.o: $(UTILDIR)/utils.a
ideq.o: $(UTILDIR)/utils.a

MODULES= $(UTILDIR)

mods: 
	for dir in $(MODULES); do (cd $$dir; $(MAKE)); done

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

$(UTILDIR)/utils.a: mods
