UTILS=utils
GEO=geo

#################################################################### OVERVIEW
#
#  Makefile for the GS2 Gyrokinetic Stability code 
#  (requires GNU's gmake)
#
#  Makefile written by Bill Dorland, Greg Hammett and Darin Ernst
#
#  LAST UPDATE: 7/29/02
#
############################################################# WALL CLOCK TIMES
#
# WALL CLOCK TIMES for s1.in (8 modes)
#
# proc   Arc 	Seaborg
#  8    14m51s 	13m46s
# 16     8m26s	 7m00s
# 32	 ...	 3m24s
#
######################################################################## TODO
#
#  First do "module load netcdf " (on Crays) or otherwise set up your
#  compilation/library environment.
#
#  Then type:
#
#    "make" to build an optimized version of GS2
#    "make debug=on" to build debug version
#    "make static=on" to build semi-portable version
#
# This uses the shell script ./test_os to figure out which platform 
# we are running on to set the right compiler switches and library 
# paths. Or type "make CPU=C90", or "setenv CPU C90" to compile for 
# Cray C90/J90. 
#
# Platform-dependent logical links are done automatically, but to override
# them (to use shmem on the t3e for example), do:
# 
#    "make t3e_shmem"  (other link options listed at end of Makefile)
#
# Alternatively, the "debug" and/or "CPU" make options can be set as 
# environment variables, and then run make.  For example:
#
#    setenv CPU T3E ; setenv debug on ; make
#
# Or set defaults by uncommenting one of the following lines:
# 
# CPU = T3E
# debug  = on
#
# To use the (lower performance) fftw library on the T3E uncomment this line:
# FFT=fftw
#
# NEW FEATURES:
#
#   "make test_make" - print values of variables
#   "make tar" - prepare a dated tar.gz file of src directory 
#   "make distclean" - does "make clean" + removes platform links & executables
# 
#
#
# CPU OPTIONS: 
#
# LINUX    	(for MIT linux cluster with NAG f95)
# LINUX_fuj 	(for linux with fujitsu f90)
# LINUX_alpha 	(for linux with Compaq Alpha compiler)
# RS6000	(for IBM SP2 with xlf)
# T3E		(for T3E with f90)
# C90		(for C90/J90 with f90) 
# ALPHA		(for DEC ALPHA with DEC f95)
# ALPHA_NAG	(for DEC ALPHA with NAG f95)
#
####################################################################### CPU
#
# environment variable CPU is used to set platform-dependent make variables.
#
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
  F90FLAGS = -I$$MPTDIR/include -M1110,7212 -p$(UTILS) -p$(GEO)
  FLIBS = -lmpi $$NETCDF -Wl"-D permok=yes" $(UTILS)/mdslib.a 
ifeq ($(FFT),fftw)
  FLIBS += -L../fftw/lib -lfftw -lrfftw
endif
  ifneq ($(debug),on)
    F90FLAGS += -O vector3 -O aggress 
  else
    F90FLAGS += -g -R abcs -e i
    FLIBS  += -Wl"-D preset=inf" # -lmalloc 
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
  F90FLAGS =  -I/usr/local/include -M1110
  FLIBS = $$NETCDF
  ifneq ($(debug),on)
    F90FLAGS += -O vector3  
  else
    F90FLAGS += -g -R abcs -e i
    FLIBS += -Wl"-D preset=inf" # -lmalloc
  endif
endif

# J90 options so that variables are private to each process (needed for MPI):
# (I can't get this to work.  So for now drop the "-a taskcommon" switch,
# which restricts us to 1 processor on the J90.)
#  F90FLAGS =  -I/usr/local/include -M1110 -O vector3 -a taskcommon
#  FLIBS_itg = -L/usr/local/lib -lnag $$NETCDF
#
# debug options:
#  F90FLAGS =     -I/usr/local/include -M1110 -g -R abcs -e i
#  FLIBS_itg = -L/usr/local/lib -lnag $$NETCDF -lmalloc -Wl"-D preset=inf"
#                 # -lapp

# SGI Origin-2000 options:
ifeq ($(CPU),SGI)
  FC = f90
  FLIBS  =    -L/usr/pppl/lib -lnetcdf -lmpi -lscs
  PLATFORM_LINKS = origin
  F90FLAGS =  -I/usr/pppl/include -64 -r8 -TARG:platform=ip27
  ifneq ($(debug),on)
    F90FLAGS += -O -TARG:platform=ip27

# Other options tried (no more than 15% speedup):
#    F90FLAGS =  -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm

  else
    F90FLAGS += -g -DEBUG:div_check=3:trap_uninitialized=on

#       -DEBUG:div_check=3:subscript_check=on:trap_uninitialized=on
# should be the full debug options, but then the compiler bombs on some
# of the routines (even the new beta version 7.3 of the compiler).

  endif
endif

# NERSC IBM options:
ifeq ($(CPU),RS6000)
  FC = mpxlf90_r
  FC5 = h5pfc
  PLATFORM_LINKS = ibm
  F90FLAGS = -qautodbl=dbl4 -qsuffix=f=f90 -I $(UTILS) -I $(GEO) \
	-I $(NETCDF_DIR)/include 
  FLIBS = $$NETCDF $$FFTW -lfftw -lrfftw
  FLIBS5 = $$NETCDF $$FFTW -lfftw -lrfftw \
	-I/usr/common/usg/hdf5/64/1.6.4/parallel/lib \
	-L/usr/common/usg/zlib/default/lib -L/usr/common/usg/szip/default/lib \
	/usr/common/usg/hdf5/64/1.6.4/parallel/lib/libhdf5_fortran.a \
	/usr/common/usg/hdf5/64/1.6.4/parallel/lib/libhdf5.a -lsz -lz -lm
  ifneq ($(debug),on)
#    F90FLAGS += -O4
    F90FLAGS += -O3 -qarch=auto -qtune=auto
  else
    F90FLAGS += -g -qflttrap=nanq
    FLIBS    += # $$TRACE_MPIF
  endif

endif

# Dawson G5 cluster:
ifeq ($(CPU),Dawson)
  FC = xlf95
  PLATFORM_LINKS = ibm
  F90FLAGS = -qmoddir=/tmp/mabarnes -I/tmp/mabarnes \
	-qautodbl=dbl4 -qsuffix=f=f90 -I $(UTILS) -I $(GEO) \
	-I /u/local/apps/netcdf/include -I/u/local/mpi/mpilam/include 
  FLIBS = -L/u/local/apps/netcdf/lib -lnetcdf \
	-L/u/home2/nfl/FFTW/lib -lfftw -lrfftw \
	-L/u/local/mpi/mpilam/lib -llammpio -llamf77mpi -lmpi -llam
#-L/u/local/apps/fftw/lib -lfftw -lrfftw 
  ifneq ($(debug),on)
#    F90FLAGS += -O4
    F90FLAGS += -O3 -qarch=g5 -qtune=g5
  else
    F90FLAGS += -g 
    FLIBS    += # $$TRACE_MPIF
  endif

endif

# DEC alpha options:
ifeq ($(CPU),ALPHA)
  FC = f95
  FLIBS = -L/usr/local/lib -lnetcdf \
	-L/u/hammett/local/alpha/lib -lfftw -lrfftw \
	$(UTILS)/mdslib.a
#	-non_shared -L/usr/local/mdsplus_new/lib -lMdsLib -lMdsShr
# (use -non_shared because it can't map MdsLib.so at run time for some 
# reason...)
  PLATFORM_LINKS = alpha
  F90FLAGS = -I/usr/local/include -r8 -I$(UTILS) -I$(GEO)

  ifeq ($(debug),on)
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -check output_conversion 
  else
     F90FLAGS += -O -fast -w
  endif

endif

# options for NAG f95 on a DEC alpha
ifeq ($(CPU),ALPHA_NAG)
  FC = /usr/local/bin/f95
  FLIBS = -L/usr/local/lib -lnetcdf
  PLATFORM_LINKS = alpha_nag
  ifeq ($(debug),on)
    F90FLAGS = -C -hpf  -I/usr/local/include -r8 -g90 -dusty
  else
    F90FLAGS =  -C -hpf -I/usr/local/include -r8 -O -dusty
  endif
endif

# options for Linux with Fujitsu f90
ifeq ($(CPU),LINUX_fuj)
  FC = mpif90
  FLIBS = -L/usr/local/lib -lnetcdf \
	-L/u/mabarnes/fftw/lib -lfftw -lrfftw \
	$(UTILS)/mdslib.a
  PLATFORM_LINKS = linux_fuj
  F90FLAGS = -A m -C cdRR8 -I/usr/local/include -X9 \
	-static-flib -Kfast -I$(UTILS) -I$(GEO)

  ifeq ($(debug),on) 
    F90FLAGS += -g -H easu
  else
    F90FLAGS += -O -f 2004,2006,2008 -Wa,--no-warn
  endif

endif

# options for Linux with Absoft f90
ifeq ($(CPU),LINUX_abs)
  FC = mpif90
  PLATFORM_LINKS = linux_abs

#  FLIBS = $(UTILS)/mdslib.a  -L/usr/local/lib -lrfftw -lfftw 
  FLIBS = $(UTILS)/mdslib.a  -L~kthcmr/fftw -lrfftw -lfftw 
  FLIBS +=  -L/usr/lib -lnetcdf

  F90FLAGS_base = -Rp -I/usr/include -N113 -p$(UTILS) -p$(GEO) 
  F90FLAGS_0    = $(F90FLAGS_base) -YEXT_SFX= 
  F90FLAGS_1    = $(F90FLAGS_base) -YEXT_SFX=_
  F90FLAGS_2    = $(F90FLAGS_base) -YEXT_SFX=__

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

# options for Linux with Lahey lf95
ifeq ($(CPU),LINUX_lf95)
  FC = mpif90
  FLIBS = -L/usr/local/lib -lnetcdf \
	-L/usr/local/lib -lfftw -lrfftw $(UTILS)/mdslib.a
  PLATFORM_LINKS = linux_lf95
  F90FLAGS = --dbl --ml cdecl -I$(UTILS) -I$(GEO) -I/usr/local/include

  ifeq ($(static),on)
    F90FLAGS += --staticlink
  endif

  ifeq ($(debug),on) 
# -mpitrace -mpianim -mpilog
    F90FLAGS += -g --chk aesu 
#    F90FLAGS += -O -mpitrace 
  else
    F90FLAGS += -O 
  endif

endif

# options for Linux with NAG f95:
ifeq ($(CPU),LINUX)
  FC = mpif90
  FLIBS = $(UTILS)/mdslib.a \
	-L/usr/local/lib -lnetcdf \
	-L/usr/local/lib -lfftw -lrfftw
  PLATFORM_LINKS = linux
  F90FLAGS = -w -r8 -mismatch -I/opt/mpich.gccf95-2/include \
	-I/usr/local/include -I $(GEO) -I $(UTILS)
  F90FLAGS_2  = -w -mismatch -I/opt/mpich.gccf95-2/include \
	-I/usr/local/include -I $(GEO) -I $(UTILS)
  F90FLAGS_3  = -r8 -mismatch -I/opt/mpich.gccf95-2/include \
	-I/usr/local/include -I $(GEO) -I $(UTILS)

  ifeq ($(debug),on)
    F90FLAGS += -C=array -C=bits -C=dangling -C=do -C=present -C=pointer -gline 
    F90FLAGS_2 += -C=array -C=bits -C=dangling -C=do -C=present -C=pointer -gline
    F90FLAGS_3 += -C=array -C=bits -C=dangling -C=do -C=present -C=pointer -gline
  else
    F90FLAGS += -O4 
  endif

endif

# 
# SUN or DEC alpha switches?:
# -eA -g -dalign -I/usr/local/include 
# -ef #-Rabc -m0

# options for Linux on a DEC alpha with DEC/Compaq F90:
# f90 is a link to "fort", for man pages do "man fort"
# (compiler options slightly different than for Compaq F90 for Ultrix):
ifeq ($(CPU),LINUX_alpha)
  FC = f90
  FLIBS = -L/u/hammett/local/alinux/lib -lfftw -lrfftw \
	$(UTILS)/mdslib.a $(UTILS)/netcdf_stub.a
#	-L/usr/local/lib -lnetcdf
#	-L/usr/local/mdsplus/lib -lMdsLib -lMdsShr
  PLATFORM_LINKS = linux_alpha
  F90FLAGS = -I/usr/local/include -r8 -I$(UTILS) -I$(GEO)

  ifeq ($(debug),on)
     F90FLAGS += -g -assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source \
	-align dcommons -align sequence
  else
     F90FLAGS +=  -O -fast -w
  endif

endif


########################################################## MODULE DECLARATIONS

# Cray compiler automatically searches for module *.o files in present
# directory, but NAG compiler doesn't, so have to explicitly list them:

GS2MOD= constants.o prof.o mp.o gs2_layouts.o command_line.o gs2_save.o \
	text_options.o file_utils.o ran.o redistribute.o antenna_data.o antenna.o \
	gs2_reinit.o gs2_time.o convert.o fft_work.o shmem.o \
	theta_grid.o kt_grids.o dist_fn_arrays.o species.o gs2_dist_io.o \
	fields_arrays.o le_grids.o collisions.o gs2_transforms.o \
	nonlinear_terms.o fields_explicit.o \
	fields.o fields_implicit.o fields_test.o init_g.o check.o \
	dist_fn.o hyper.o gs2_diagnostics.o gs2_io.o netcdf_mod.o \
	run_parameters.o gs2_flux.o regression.o gs2_heating.o \
	$(GEO)/geo.a gridgen4mod.o $(UTILS)/spl.o $(UTILS)/utils.a
#	utils/utils.a geo/geo.a utils/mds.o 
# *.a libraries must appear after any function that needs them.

INGENMOD= constants.o mp.o gs2_layouts.o command_line.o \
	text_options.o file_utils.o theta_grid.o \
	$(GEO)/geo.a gridgen4mod.o $(UTILS)/spl.o $(UTILS)/utils.a

EGRIDMOD = $(UTILS)/spl.o constants.o

ifeq ($(CPU),LINUX)
  GS2MOD += nag_args.o
  mp.o: nag_args.o
endif

# ifneq ($(loc),mit) 
#	GS2MOD += utils/mdslib.o 
# endif

LINKS= command_line.f90 mp.f90 shmem.f90 prof.f90 redistribute.f90 ran.f90 gs2_layouts.f90 \
	gs2_save.f90 gs2_transforms.f90 fft_work.f90 $(UTILS)/mds.f90 check.f90

######################################################################## RULES
.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o: 
	$(FC) $(F90FLAGS) -c $<

################################################################### DIRECTIVES
# search path where gnu make will look for files:
# VPATH = ./src

# check links and subdirectory modules before building gs2 itself:
# (This order is important because some dependencies are known only
# by the Makefiles in the subdirectories, etc.).
all: $(LINKS) modules gs2 ingen rungridgen

gs2: gs2.o $(GS2MOD) 
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(F90FLAGS) -o gs2 gs2.o $(UTILS)/mpptime.o $(FLIBS) ;; \
		ibm) $(FC) $(F90FLAGS) -o gs2 gs2.o $(GS2MOD) $(FLIBS5) ;; \
		*)   $(FC) $(F90FLAGS) -o gs2 gs2.o $(GS2MOD) $(FLIBS) ;; \
	esac

gs2.x: gs2.o $(GS2MOD) 
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(F90FLAGS) -o gs2.x gs2.o $(UTILS)/mpptime.o $(FLIBS) ;;\
		*)   $(FC) $(F90FLAGS) -o gs2.x gs2.o $(GS2MOD) $(FLIBS) ;; \
	esac

ingen: ingen.o $(INGENMOD)
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(F90FLAGS) -o ingen ingen.o $(FLIBS) ;; \
		*)   $(FC) $(F90FLAGS) -o ingen ingen.o $(INGENMOD) $(FLIBS) ;; \
	esac

egrid: egrid.o $(EGRIDMOD)
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(F90FLAGS) -o egrid egrid.o $(FLIBS) ;; \
		*)   $(FC) $(F90FLAGS) -o egrid egrid.o $(EGRIDMOD) $(FLIBS) ;; \
	esac

regress: regression.o drive.o
	$(FC) $(F90FLAGS) -o regress drive.o regression.o $(FLIBS)

rgobj = gridgen4mod.o file_utils.o mp.o command_line.o $(UTILS)/utils.a 

rungridgen: rungridgen.o gridgen4mod.o file_utils.o command_line.o $(UTILS)/utils.a
	case $(PLATFORM_LINKS) in \
		t3e) $(FC) $(F90FLAGS) -o rungridgen rungridgen.o ;; \
		*)   $(FC) $(F90FLAGS) -o rungridgen rungridgen.o $(rgobj) ;; \
	esac

# libraries and functions in subdirectories:
# ??? Code bombs on Ultrix Alpha on test case s1a.in...

modules: 
	@cd $(UTILS) ; $(MAKE)
	@cd $(GEO)   ; $(MAKE)

$(UTILS)/utils.a:
	cd $(UTILS) ; $(MAKE)

$(GEO)/geo.a:
	cd $(GEO) ; $(MAKE)

#utils/spl.o: 
#	cd utils; $(MAKE)

file_utils.o: 
	$(FC) $(F90FLAGS) -c file_utils.f90

gs2_dist_io.o: gs2_dist_io.f90
	$(FC5) $(F90FLAGS) -c gs2_dist_io.f90

ifeq ($(CPU),LINUX_lf95)

le_grids.o: le_grids.f90
	$(FC) --quad $(F90FLAGS) -c le_grids.f90
endif

ifeq ($(CPU),Dawson)

le_grids.o: le_grids.f90
	$(FC) $(F90FLAGS) -qautodbl=dbl -c le_grids.f90
endif

ifeq ($(CPU),RS6000)

le_grids.o: le_grids.f90
	$(FC) $(F90FLAGS) -qautodbl=dbl -c le_grids.f90
endif

#
# No optimizations for some routines: stupid workarounds because of 
# incompatibility between f95 and my C libraries (I think)
ifeq ($(CPU),LINUX)

gs2_transforms.o: gs2_transforms.f90
	$(FC) $(F90FLAGS_3) -c $<

nonlinear_terms.o: nonlinear_terms.f90
	$(FC) $(F90FLAGS_3) -c $<

netcdf_mod.o: netcdf_mod.f90
	$(FC) $(F90FLAGS_2) -c netcdf_mod.f90

endif
ifeq ($(CPU),LINUX_abs)

gs2_io.o:
	$(FC) $(F90FLAGS_2) -c gs2_io.f90
gs2_save.o:
	$(FC) $(F90FLAGS_2) -c gs2_save.f90
netcdf_mod.o:
	$(FC) $(F90FLAGS_2) -c netcdf_mod.f90

command_line.o:
	$(FC) $(F90FLAGS_1) -c command_line.f90
file_utils.o: 
	$(FC) $(F90FLAGS_1) -c file_utils.f90
mp.o:
	$(FC) $(F90FLAGS_1) -c mp.f90

fft_work.o:
	$(FC) $(F90FLAGS_0) -c fft_work.f90 
gs2_transforms.o:  
	$(FC) $(F90FLAGS_0) -c gs2_transforms.f90

endif

################################################################# DEPENDENCIES

file_utils.o: command_line.o mp.o file_utils.f90
antenna.o: species.o run_parameters.o file_utils.o mp.o gs2_time.o antenna_data.o
antenna.o: kt_grids.o theta_grid.o ran.o constants.o 
gs2_layouts.o: mp.o file_utils.o
run_parameters.o: mp.o file_utils.o gs2_save.o kt_grids.o text_options.o 
run_parameters.o: gs2_time.o
species.o: mp.o file_utils.o text_options.o
gs2_save.o: theta_grid.o gs2_layouts.o mp.o fields_arrays.o kt_grids.o 
gs2_save.o: file_utils.o antenna_data.o
gs2_transforms.o: gs2_layouts.o mp.o prof.o fft_work.o redistribute.o
gs2_transform.o: theta_grid.o kt_grids.o 
gs2_diagnostics.o: file_utils.o kt_grids.o run_parameters.o species.o mp.o 
gs2_diagnostics.o: fields.o dist_fn.o constants.o prof.o gs2_save.o gs2_time.o
gs2_diagnostics.o: gs2_io.o le_grids.o fields_arrays.o dist_fn_arrays.o 
gs2_diagnostics.o: gs2_transforms.o nonlinear_terms.o collisions.o $(UTILS)/utils.a
gs2_diagnostics.o: gs2_flux.o gs2_heating.o gs2_dist_io.o
gs2_heating.o: mp.o species.o
dist_fn.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o antenna.o
dist_fn.o: run_parameters.o init_g.o text_options.o fft_work.o gs2_heating.o
dist_fn.o: gs2_layouts.o file_utils.o dist_fn_arrays.o constants.o  gs2_time.o
dist_fn.o: collisions.o nonlinear_terms.o
dist_fn.o: gs2_transforms.o prof.o gs2_time.o redistribute.o hyper.o
hyper.o: kt_grids.o run_parameters.o file_utils.o text_options.o 
hyper.o: mp.o gs2_layouts.o theta_grid.o gs2_time.o le_grids.o
init_g.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o dist_fn_arrays.o
init_g.o: gs2_layouts.o gs2_save.o fields_arrays.o ran.o text_options.o
init_g.o: file_utils.o run_parameters.o
le_grids.o: mp.o species.o theta_grid.o kt_grids.o file_utils.o redistribute.o
le_grids.o: gs2_layouts.o constants.o 
gs2.o: mp.o file_utils.o fields.o run_parameters.o gs2_diagnostics.o 
gs2.o: gs2_reinit.o gs2_time.o init_g.o check.o dist_fn_arrays.o
gs2.o: gs2_save.o 
fields.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields.o: file_utils.o dist_fn_arrays.o constants.o prof.o text_options.o
fields.o: fields_arrays.o fields_implicit.o fields_test.o fields_explicit.o 
fields.o: init_g.o nonlinear_terms.o 
fields_implicit.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields_implicit.o: fields_arrays.o dist_fn_arrays.o antenna.o
fields_implicit.o: gs2_layouts.o prof.o file_utils.o
fields_implicit.o: gs2_save.o nonlinear_terms.o
fields_explicit.o: theta_grid.o kt_grids.o dist_fn.o mp.o
fields_explicit.o: fields_arrays.o dist_fn_arrays.o gs2_save.o 
fields_explicit.o: run_parameters.o gs2_layouts.o prof.o file_utils.o
fields_test.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields_test.o: fields_arrays.o dist_fn_arrays.o
fields_test.o: gs2_layouts.o mp.o prof.o file_utils.o
kt_grids.o: mp.o file_utils.o text_options.o theta_grid.o constants.o 
theta_grid.o: mp.o file_utils.o text_options.o constants.o gridgen4mod.o $(GEO)/geo.a
theta_grid.o: $(UTILS)/utils.a
collisions.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o gs2_time.o
collisions.o: run_parameters.o file_utils.o text_options.o dist_fn_arrays.o
collisions.o: prof.o shmem.o redistribute.o gs2_layouts.o constants.o 
nonlinear_terms.o: theta_grid.o kt_grids.o le_grids.o species.o gs2_layouts.o 
nonlinear_terms.o: dist_fn_arrays.o gs2_transforms.o run_parameters.o constants.o 
nonlinear_terms.o: text_options.o mp.o gs2_time.o file_utils.o 
gs2_reinit.o: collisions.o mp.o nonlinear_terms.o gs2_time.o
gs2_reinit.o: fields_explicit.o dist_fn.o fields.o fields_implicit.o fields_test.o
gs2_reinit.o: init_g.o run_parameters.o gs2_save.o dist_fn_arrays.o fields_arrays.o
gs2_reinit.o: file_utils.o antenna.o  #additional_terms.o
redistribute.o: mp.o
gs2_io.o: mp.o file_utils.o netcdf_mod.o kt_grids.o theta_grid.o le_grids.o 
gs2_io.o: species.o run_parameters.o convert.o fields_arrays.o 
gs2_io.o: nonlinear_terms.o gs2_layouts.o constants.o gs2_transforms.o gs2_heating.o
check.o: mp.o file_utils.o run_parameters.o 
netcdf_mod.o: mp.o constants.o
ingen.o: file_utils.o text_options.o constants.o theta_grid.o
gs2_flux.o: species.o mp.o text_options.o file_utils.o dist_fn.o regression.o 
gridgen4mod.o:  $(UTILS)/utils.a
gs2_dist_io.o: mp.o gs2_transforms.o kt_grids.o gs2_layouts.o theta_grid.o
gs2_dist_io.o: le_grids.o species.o file_utils.o

############################################################## MORE DIRECTIVES
clean:
	rm -f *.o *.mod *.g90 core */core ; \
	cd $(UTILS) ; rm -f *.o *.mod *.a ; cd .. ; \
	cd $(GEO) ; rm -f *.o *.mod *.a ; cd .. ; \

distclean:	unlink clean
	rm -f rungridgen
	rm -f $(GEO)/ball
	rm -f $(GEO)/eiktest
	rm -f gs2
	rm -f ingen

tar: 
	cd ..; \
        ls src_$(DATE) || ln -s src src_$(DATE); \
        tar cvf GS2_src_$(DATE).tar \
           src_$(DATE)/*.f* \
           src_$(DATE)/Inputs/*.in \
           src_$(DATE)/Outputs/*.out \
           src_$(DATE)/*.in \
           src_$(DATE)/README.ernst \
           src_$(DATE)/Makefile* \
           src_$(DATE)/*.inf \
           src_$(DATE)/*.html  \
           src_$(DATE)/test*os ; \
        gzip GS2_src_$(DATE).tar 


############################################################### DEFAULT RULES

# If no other rules are found, use the defaults:

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

test_make:
	@echo FC is $(FC)
	@echo F90FLAGS is $(F90FLAGS)
	@echo FLIBS is $(FLIBS)
	@echo debug is $(debug)
	@echo CPU is $(CPU)

# If one of the platform-specific logical links doesn't exist, set them up:

$(LINKS):
	@$(MAKE) --no-print-directory $(PLATFORM_LINKS)

############################################################### PLATFORM LINKS
#
# Platform-specific logical links:
#

c90:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_sgi.f90 check.f90
	ln -sf ran_cray.f90 ran.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_unicos.f90 fft_work.f90
	ln -sf gs2_transforms_sgi.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

t3e_shmem:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi.f90 mp.f90
	ln -sf shmem_cray.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_shnew.f90 redistribute.f90
	ln -sf check_sgi.f90 check.f90
	ln -sf ran_cray.f90 ran.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_unicosmk.f90 fft_work.f90
	ln -sf gs2_transforms_sgi.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

t3e:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_sgi.f90 check.f90
	ln -sf ran_cray.f90 ran.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_unicosmk.f90 fft_work.f90
	ln -sf gs2_transforms_sgi.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

t3e_fftw:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_shnew.f90 redistribute.f90
	ln -sf check_sgi.f90 check.f90
	ln -sf ran_cray.f90 ran.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

ibm:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_aix.f90 check.f90
	ln -sf ran_portable.f90 ran.f90
	ln -sf file_utils_xlf.f90 file_utils.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_save_aix.f90 gs2_save.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

origin:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_sgi.f90 check.f90
	ln -sf ran_cray.f90 ran.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_origin.f90 fft_work.f90
	ln -sf gs2_transforms_sgi.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

linux:
	ln -sf command_line_nag.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_portable.f90 ran.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd $(UTILS); ln -sf mds_io_stub.f90 mds.f90 

linux_fuj:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_local.f90 ran.f90
#	ln -sf gs2_save_stub.f90 gs2_save.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..
	ln -sf gs2_transforms_stub.f90 gs2_transforms.f90
	ln -sf fft_work_stub.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90

linux_lf95:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_portable.f90 ran.f90
#	ln -sf gs2_save_stub.f90 gs2_save.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90

linux_abs:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_local.f90 ran.f90
#	ln -sf gs2_save_stub.f90 gs2_save.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90

alpha:
	ln -sf command_line_alpha.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf ran_local.f90 ran.f90
	ln -sf check_portable.f90 check.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..
#	cd utils; ln -sf mds_io.f90 mds.f90 ; cd ..

alpha_nag:
	ln -sf command_line_nag.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_local.f90 ran.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..

linux_alpha:
	ln -sf command_line_alpha.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_local.f90 ran.f90
#	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf gs2_save_stub.f90 gs2_save.f90
#	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_transforms_stub.f90 gs2_transforms.f90
	ln -sf file_utils_portable.f90 file_utils.f90
#	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf fft_work_stub.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	cd utils; ln -sf mds_io_stub.f90 mds.f90 ; cd ..

unlink:
	rm -f command_line.f90
	rm -f mp.f90
	rm -f shmem.f90
	rm -f prof.f90
	rm -f redistribute.f90
	rm -f ran.f90
	rm -f gs2_save.f90
	rm -f gs2_transforms.f90
	rm -f fft_work.f90
	rm -f $(UTILS)/mds.f90
	rm -f check.f90 
	rm -f file_utils.f90
	rm -f gs2_dist_io.f90

