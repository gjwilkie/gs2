#################################################################### OVERVIEW
#
#  Makefile for the GS2/AstroGK Gyrokinetic Turbulence code 
#  (requires GNU's gmake)
#
PROJECT ?= gs2
#PROJECT ?= agk
#
#  Makefile written by Bill Dorland and Ryusuke Numata
#
#  LAST UPDATE: 06/16/08
#
# * Changelogs
#	07/01/08: new switches USE_LOCAL_SPFUNC and USE_NAGLIB
#	06/16/08: clean up unused statements related to PLATFORM_LINKS
#
# * Available Compilers (tested on limited hosts)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# --- NAGWare
# --- Fujitsu
# --- Absoft
# --- Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (... with Intel)
#
# * Switches:
#
# Here is a list of switches with simple explanation.
# In the brackets, values accepted are shown,
# where "undefined" means blank.
# Switches with (bin) are checked if they are defined or not
# What values they have do not matter.
# Be careful that DEBUG=off means DEBUG=on.
#
# turns on debug mode (bin)
DEBUG ?=
# turns on test mode (bin)
TEST ?=
# turns on profile mode (gprof,ipm)
PROF ?=
# optimization (on,aggressive,undefined)
OPT ?= on
# prevents linking with shared libraries (bin)
STATIC ?=
# promotes precisions of real and complex (bin)
DBLE ?= on
# turns on distributed memory parallelization using MPI (bin)
USE_MPI ?= on
# turns on SHMEM parallel communications on SGI (bin)
USE_SHMEM ?=
# which FFT library to use (fftw,undefined) 
USE_FFT ?= fftw
# uses netcdf library (new,old,undefined)
USE_NETCDF ?= new
# uses hdf5 library (bin)
USE_HDF5 ?=
# uses MDSplus (bin)
USE_MDSPLUS ?=
# Use function pointer in agk_layouts_indices.c (bin)
# see also README.cpp
USE_C_INDEX ?= 
# Use Numerical Recipes local random number generator (bin)
# see also README.cpp
USE_NR_RAN ?=
# Use posix for command_line (bin)
USE_POSIX ?=
# Use local special functions (bin)
USE_LOCAL_SPFUNC ?= 
# Use nag libraray (spfunc,undefined)
USE_NAGLIB ?= 
#
# * Targets:
#
#  depend: generate dependency
#  test_make: print values of variables
#  clean: clean up
#  distclean: does "make clean" + removes platform links & executables
#  tar: pack
#
############################################################### DEFAULT VALUE
#
# These variables can be set in platform-dependent Makefile.
#

MAKE = make
CPP = cpp
CPPFLAGS = -C -P -traditional
FC = f90
F90FLAGS =
F90OPTFLAGS =
CC = cc
CFLAGS =
COPTFLAGS =
LD = $(FC)
LDFLAGS = $(F90FLAGS)
ARCH = ar
ARCHFLAGS = cr
RANLIB=ranlib
AWK = awk

# These macros are used for the suffix problem of absoft
F90FLAGS_SFX0 =
F90FLAGS_SFX1 =
F90FLAGS_SFX2 =

MPI_INC =
MPI_LIB =
FFT_INC =
FFT_LIB =
NETCDF_INC =
NETCDF_LIB =
HDF5_INC =
HDF5_LIB =
IPM_LIB =
NAG_LIB =
NAG_PREC ?= dble
PGPLOT_LIB =

################################################### SET COMPILE MODE SWITCHES

ifdef TEST
$(warning TEST mode is not working yet)
	override TEST =
endif

######################################################### PLATFORM DEPENDENCE

# compile mode switches (DEBUG, TEST, PROF, OPT, STATIC, DBLE)
# must be set before loading Makefile.$(SYSTEM) because they may affect
# compiler options.
# However, Makefile.local may override some options set in Makefile.$(SYSTEM),
# thus it is included before and after Makefile.$(SYSTEM)
sinclude Makefile.local

# include system-dependent make variables
sinclude Makefile.$(SYSTEM)

# include Makefile.local if exists
sinclude Makefile.local

#############################################################################

UTILS=utils
GEO=geo

ifeq ($(MAKECMDGOALS),depend)
# must invoke full functionality when make depend
	MAKE += USE_HDF5=on USE_FFT=fftw USE_NETCDF=new USE_MPI=on \
		USE_LOCAL_BESSEL=on
endif

ifdef USE_SHMEM
$(warning USE_SHMEM is not working yet)
	override USE_SHMEM =	
endif

ifdef USE_HDF5
	ifndef USE_MPI
$(error Currently, USE_HDF5 works with USE_MPI)
	endif
endif

ifndef USE_FFT
$(warning USE_FFT is off)
$(warning Be sure that nonlinear run makes no sense)
endif

ifdef USE_MPI
	CPPFLAGS += -DMPI
endif
ifdef USE_SHMEM
	CPPFLAGS += -DSHMEM
endif
ifeq ($(USE_FFT),fftw)
	CPPFLAGS += -DFFT=_FFTW_
	FFT_LIB ?= -lfftw -lrfftw
endif
ifdef USE_NETCDF
	ifeq ($(USE_NETCDF),new)
		CPPFLAGS += -DNETCDF=_DEFAULT_
	endif
	ifeq ($(USE_NETCDF),old)
		CPPFLAGS += -DNETCDF=_OLD_
	endif
	NETCDF_LIB ?= -lnetcdf
endif
ifdef USE_HDF5
	ifdef USE_MPI
		FC = h5pfc
	endif
	CPPFLAGS += -DHDF
endif
ifdef USE_C_INDEX
	CPPFLAGS += -DUSE_C_INDEX
endif
ifdef USE_NR_RAN
	CPPFLAGS += -DUSE_NR_RAN
endif
ifdef USE_POSIX
	CPPFLAGS += -DPOSIX
endif
ifdef USE_LOCAL_SPFUNC
	CPPFLAGS += -DSPFUNC=_SPLOCAL_
else
	ifeq ($(findstring spfunc,$(USE_NAGLIB)),spfunc)
		CPPFLAGS += -DSPFUNC=_SPNAG_
	endif
endif
ifdef USE_NAGLIB
	ifeq ($(NAG_PREC),dble)
		ifndef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGDBLE_
	endif
	ifeq ($(NAG_PREC),sngl)
		ifdef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		CPPFLAGS += -DNAG_PREC=_NAGSNGL_
	endif
endif
ifndef PGPLOT_LIB
	ifeq ($(MAKECMDGOALS),agk_fields_plot)
$(error PGPLOT_LIB is not defined)
	endif
endif

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(IPM_LIB) $(NAG_LIB)
PLIBS 	+= $(LIBS) $(PGPLOT_LIB)
F90FLAGS+= $(F90OPTFLAGS) \
	   $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC)
CFLAGS += $(COPTFLAGS)

DATE=$(shell date +%y%m%d)
TARDIR=$(PROJECT)_$(DATE)
TOPDIR=$(CURDIR)
ifeq ($(notdir $(CURDIR)), $(UTILS))
	TOPDIR=$(subst /$(UTILS),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(GEO))
	TOPDIR=$(subst /$(GEO),,$(CURDIR))
endif
ifneq ($(TOPDIR),$(CURDIR))
	SUBDIR=true
endif

VPATH = $(UTILS):$(GEO):Aux:../$(UTILS):../$(GEO)
# this just removes non-existing directory from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH = .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")
#
ifdef SUBDIR
	VPATH +=:..
endif
DEPEND=Makefile.depend
PERL_CMD=perl
DEPEND_CMD=$(PERL_CMD) fortdep

# most common include and library directories
DEFAULT_INC_LIST = . $(UTILS) $(GEO) .. ../$(UTILS) ../$(GEO) \
		/usr/include /usr/local/include \
	   	/opt/local/include /sw/include
DEFAULT_LIB_LIST =
#DEFAULT_LIB_LIST = /usr/lib /usr/local/lib \
#		/opt/local/lib /sw/lib
# This default library path list can simplify the procedure of porting,
# however, I found this (actually -L/usr/lib flag) causes an error
# when linking gs2 at bassi (RS6000 with xl fortran)
DEFAULT_INC=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

####################################################################### RULES

.SUFFIXES:
.SUFFIXES: .fpp .f90 .c .o

.f90.o: 
	$(FC) $(F90FLAGS) -c $<
.fpp.f90:
	$(CPP) $(CPPFLAGS) $< $@
.c.o:
	$(CC) $(CFLAGS) -c $<

##################################################################### TARGETS

# .DEFAULT_GOAL works for GNU make 3.81 (or higher)
# For 3.80 or less, see all target
.DEFAULT_GOAL := $(PROJECT)_all
ifeq ($(notdir $(CURDIR)),utils)
	.DEFAULT_GOAL := utils_all
endif
ifeq ($(notdir $(CURDIR)),geo)
	.DEFAULT_GOAL := geo_all
endif

.PHONY: all gs2_all agk_all

all: $(.DEFAULT_GOAL)

include $(DEPEND)

gs2_all: modules gs2 ingen rungridgen

gs2: $(gs2_mod) 
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

gs2.x: $(gs2_mod) 
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

agk_all: agk

ifeq ($(PROJECT),agk)
ifdef USE_C_INDEX
astrogk_mod += agk_layouts_indices.o
endif
endif

agk: $(astrogk_mod) 
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

agk.x: $(astrogk_mod) 
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

ingen: $(ingen_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

rungridgen: $(rungridgen_mod)
	$(LD) $(LDFLAGS) -o $@ $^

regress: $(drive_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

slice_g: $(slice_g_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

agk_fields_plot: $(agk_fields_plot_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(PLIBS)

.PHONY: modules utils_all geo_all

modules: utils.a geo.a

### UTILS
utils_all: utils.a mdslib.a

UTIL_OBJ = spl.o mds.o constants.o file_utils.o netcdf_utils.o
utils.a: $(UTIL_OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

mdslib.a: mdslib.o
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

### GEO
geo_all: geo.a eiktest ball

GEO_OBJ = geometry.o geq.o veq.o eeq.o peq.o leq.o deq.o ideq.o \
	radstub.o vdimstub.o mdslib.o
geo.a: $(GEO_OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

ball: $(ball_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

eiktest: $(eiktest_mod)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

############################################################### SPECIAL RULES

# These are special rules for the suffix problem of absoft
# (not tested)
$(PROJECT)_transforms.o: $(PROJECT)_transforms.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $*.f90
	rm -f $*.f90
$(PROJECT)_io.o: $(PROJECT)_io.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $*.f90
	rm -f $*.f90
$(PROJECT)_save.o: $(PROJECT)_save.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $*.f90
	rm -f $*.f90
mp.o: mp.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $*.f90
	rm -f $*.f90
fft_work.o: fft_work.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $*.f90
	rm -f $*.f90

agk_layouts_indices.o: agk_layouts_type.h
agk_layouts_type.h: agk_layouts_type.f90
	$(AWK) -f makehead.awk $^ > $@

############################################################# MORE DIRECTIVES

.PHONY: depend clean distclean tar test_make

depend:	
	@$(DEPEND_CMD) -m "$(MAKE)" -1 -o -v=0 $(VPATH)

clean:
	-rm -f *.o *.mod *.g90 *.h core */core

distclean: unlink clean
	-rm -f $(PROJECT)
	-rm -f *.a
	-rm -f ingen rungridgen regress
	-rm -f ball eiktest
	-rm -f slice_g
	-rm -f agk_fields_plot

tar:
	@[ ! -d $(TARDIR) ] || echo "ERROR: directory $(TARDIR) exists. Stop."
	@[ -d $(TARDIR) ] || $(MAKE) tar_exec

### setting tar_exec local $(TARLIST*) variables
# expand wildcards listed $(TARLIST_wild) in ( $(TARLIST_dir) + . )
# directories and add them into TARLIST
tar_exec: TARLIST = test_os makehead.awk fortdep AstroGK.in
tar_exec: TARLIST_dir = Makefiles utils geo Aux
tar_exec: TARLIST_wild = *.f90 *.fpp *.inc *.c Makefile Makefile.* README README.*
tar_exec: TARLIST += $(foreach dir,. $(TARLIST_dir),$(wildcard $(addprefix $(dir)/,$(TARLIST_wild))))

tar_exec:
	@mkdir $(TARDIR)
	@for dir in $(TARLIST_dir) ;\
	  do ( [ ! -d $$dir ] ||  mkdir $(TARDIR)/$$dir ; ) ;\
	done
	@for name in $(TARLIST) ;\
	  do ( [ -f $$name ] && ln $$name $(TARDIR)/$$name ; ) ;\
	done
	@tar cvf - $(TARDIR) | bzip2 -9 > $(TARDIR).tar.bz2
	@rm -rf $(TARDIR)

test_make:
	@echo SYSTEM is $(SYSTEM)
	@echo .DEFAULT_GOAL is $(.DEFAULT_GOAL)
	@echo VPATH is $(VPATH)
	@echo CURDIR is $(CURDIR)
	@echo TOPDIR is $(TOPDIR)
	@echo
	@echo Compile mode:
	@echo  DEBUG is $(DEBUG)
	@echo  TEST is $(TEST)
	@echo  PROF is $(PROF)
	@echo  OPT is $(OPT)
	@echo  STATIC is $(STATIC)
	@echo  DBLE is $(DBLE)
	@echo
	@echo Functions:
	@echo  USE_MPI is $(USE_MPI)
	@echo  USE_SHMEM is $(USE_SHMEM)
	@echo  USE_FFT is $(USE_FFT)
	@echo  USE_NETCDF is $(USE_NETCDF)
	@echo  USE_HDF5 is $(USE_HDF5)
	@echo  USE_MDSPLUS is $(USE_MDSPLUS)
	@echo  USE_C_INDEX is $(USE_C_INDEX)
	@echo  USE_NR_RAN is $(USE_NR_RAN)
	@echo  USE_POSIX is $(USE_POSIX)
	@echo  USE_LOCAL_SPFUNC is $(USE_LOCAL_SPFUNC)
	@echo  USE_NAGLIB is $(USE_NAGLIB)
	@echo
	@echo FC is $(FC)
	@echo F90FLAGS is $(F90FLAGS)
	@echo F90OPTFLAGS is $(F90OPTFLAGS)
	@echo CC is $(CC)
	@echo CFLAGS is $(CFLAGS)
	@echo COPTFLAGS is $(COPTFLAGS)
	@echo LD is $(LD)
	@echo LDFLAGS is $(LDFLAGS)
	@echo CPP is $(CPP)
	@echo CPPFLAGS is $(CPPFLAGS)
	@echo LIBS is $(LIBS)
	@echo PLIBS is $(PLIBS)

unlink:
	@for name in `find . -name "*.fpp"` ;\
	do \
	 namef90=`basename $$name | sed "s/\.fpp/\.f90/g"` ;\
	 if [ -f $$namef90 ]; then \
	  echo "rm $$namef90" ;\
	  rm -f $$namef90 ;\
	 fi ;\
	done
	-rm -f agk_layouts_type.h
