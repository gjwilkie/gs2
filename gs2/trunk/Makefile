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
#  LAST UPDATE: 04/15/08
#
# * TODO
#
# Write a README file which fully explains new make system.
# Replace PLATFORM_LINKS by cpp mechanism in gs2. 
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
# which FFT library to use (fftw,undefined) 
USE_FFT ?= fftw
# uses netcdf library (new,old,undefined)
USE_NETCDF ?= new
# uses hdf5 library (bin)
USE_HDF5 ?=
# Use function pointer in agk_layouts_indices.c (bin)
# see also README.cpp
USE_C_INDEX ?= 
# Use Numerical Recipes local random number generator (bin)
# see also README.cpp
USE_NR_RAN ?=
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

# work around
F90FLAGS_real_double_promote = $(F90FLAGS)
F90FLAGS_no_real_promotion = $(F90FLAGS)
F90FLAGS_without_minus_w = $(F90FLAGS)
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

PLATFORM_LINKS ?= undefined

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
	MAKE += USE_HDF5=on USE_FFT=fftw USE_NETCDF=new USE_MPI=on
endif

ifdef TEST
$(warning TEST mode is not working yet)
endif

ifdef USE_HDF5
	ifndef USE_MPI
$(error Currently, USE_HDF5 works with USE_MPI)
	endif
endif
ifeq ($(PROJECT),gs2) 
	ifndef DBLE
$(error DBLE cannot be off for gs2)
	endif
	ifndef USE_FFT
$(error USE_FFT cannot be off for gs2)
		USE_FFT=fftw
	endif
	ifndef USE_NETCDF
$(error USE_NETCDF cannot be off for gs2)
		USE_NETCDF=new
	endif
endif

ifdef USE_MPI
	CPPFLAGS += -DMPI
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

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB) \
		$(IPM_LIB)
#LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(FFT_LIB) $(NETCDF_LIB) $(HDF5_LIB)
F90FLAGS+= $(F90OPTFLAGS) \
	   $(DEFAULT_INC) $(MPI_INC) $(FFT_INC) $(NETCDF_INC) $(HDF5_INC)
CFLAGS += $(COPTFLAGS)

# TT> format changed to yymmdd
#DATE=$(shell date +%D | sed 's/\///g')
DATE=$(shell date +%y%m%d)
# <TT
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

######################################################### MODULE DECLARATIONS

ifeq ($(PROJECT),gs2)
LINKS_gs2 = \
	command_line.f90 mp.f90 shmem.f90 prof.f90 redistribute.f90 ran.f90 \
	gs2_layouts.f90 gs2_save.f90 gs2_transforms.f90 fft_work.f90 \
	check.f90
	ifdef SUBDIR
		LINKS_gs2 += ../mds.f90
	else
		LINKS_gs2 += mds.f90
	endif
endif

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

ifeq ($(PROJECT),gs2)
ifeq ($(PLATFORM_LINKS),undefined)
$(error PLATFORM_LINKS is $(PLATFORM_LINKS))
endif
endif

.DEFAULT_GOAL = $(PROJECT)_all
ifeq ($(notdir $(CURDIR)),utils)
	.DEFAULT_GOAL = utils_all
endif
ifeq ($(notdir $(CURDIR)),geo)
	.DEFAULT_GOAL = geo_all
endif

.PHONY: all gs2_all agk_all

all: $(.DEFAULT_GOAL)

include $(DEPEND)

gs2_all: $(LINKS_gs2) modules gs2 ingen rungridgen

ifeq ($(PROJECT),gs2)
# define module in mdslib.f90, otherwise the script cannot determine
# dependency correctly
gs2_mod += mdslib.o
ingen_mod += mdslib.o
eiktest_mod += mdslib.o
ball_mod += mdslib.o
eeq.o: mdslib.o
endif

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

# libraries and functions in subdirectories:
# ??? Code bombs on Ultrix Alpha on test case s1a.in...

.PHONY: modules utils_all geo_all

modules: utils.a geo.a

### UTILS
utils_all: utils.a mdslib.a

UTIL_OBJ = spl.o mds.o prec.o netcdf.o
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

ifeq ($(PROJECT),gs2)
# le_grids needs quad precision on some hosts
le_grids.o: le_grids.f90
	$(FC) $(F90FLAGS_real_double_promote) -c $<
# mds neeeds quad precision on some hosts
mds.o: mds.f90
	$(FC) $(F90FLAGS_real_double_promote) -c $<
endif
# No optimizations for some routines: stupid workarounds because of 
# incompatibility between f95 and my C libraries (I think)

# NAG: without -w
gs2_transforms.o: gs2_transforms.f90
	$(FC) $(F90FLAGS_without_minus_w) $(F90FLAGS_SFX0) -c $<
agk_transforms.o: agk_transforms.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS_without_minus_w) $(F90FLAGS_SFX0) -c $*.f90
	rm -f $*.f90
nonlinear_terms.o: nonlinear_terms.f90
	$(FC) $(F90FLAGS_without_minus_w) -c $<

# NAG: without -r8
netcdf_mod.o: netcdf_mod.f90
	$(FC) $(F90FLAGS_no_real_promotion) $(F90FLAGS_SFX2) -c $<

# Absoft: needs special treatment for suffix(?)
gs2_io.o: gs2_io.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
gs2_save.o: gs2_save.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
agk_io.o: agk_io.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $*.f90
	rm -f $*.f90
agk_save.o: agk_save.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $*.f90
	rm -f $*.f90
#netcdf_mod.o: netcdf_mod.f90
#	$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<

ifeq ($(PROJECT),gs2)
command_line.o: command_line.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $<
file_utils.o: file_utils.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $<
mp.o: mp.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $<
fft_work.o: fft_work.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<
endif
ifeq ($(PROJECT),agk)
command_line.o: command_line.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $*.f90
	rm -f $*.f90
file_utils.o: file_utils.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $*.f90
	rm -f $*.f90
mp.o: mp.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $*.f90
	rm -f $*.f90
fft_work.o: fft_work.fpp
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $*.f90
	rm -f $*.f90
endif
#gs2_transforms.o: gs2_transforms.f90
#	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<
#agk_transforms.o: agk_transforms.f90
#	$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<

agk_layouts_indices.o: agk_layouts_type.h
agk_layouts_type.h: agk_layouts_type.f90
	$(AWK) -f makehead.awk $^ > $@

############################################################# MORE DIRECTIVES

.PHONY: depend clean distclean tar test_make

# this dependency does not appear for some PLATFORM_LINKS
check.o: run_parameters.o
ifeq ($(PROJECT),gs2)
depend:
	@$(MAKE) -C $(TOPDIR) $(PLATFORM_LINKS)
	@$(DEPEND_CMD) -1 -o -v=0 $(VPATH)
else
depend:	
	@$(DEPEND_CMD) -m "$(MAKE)" -1 -o -v=0 $(VPATH)
endif

clean:
	-rm -f *.o *.mod *.g90 *.h core */core

distclean: unlink clean
	-rm -f $(PROJECT)
	-rm -f *.a
	-rm -f ingen rungridgen regress
	-rm -f ball eiktest
	-rm -f slice_g

tar: tar_$(PROJECT)

tar_gs2:
	@echo $(PROJECT)_$(DATE) > .package
	@-rm -rf `cat .package`
	@mkdir `cat .package`
	@for dir in utils geo Makefiles ;\
	do ( mkdir `cat .package`/$$dir ; ) ;\
	done
	@for name in *.f90 test_os Makefile Makefile.* \
		utils/*.f90 utils/Makefile \
		geo/*.f90 geo/Makefile \
		Makefiles/*; \
	do (ln $$name `cat .package`/$$name ; ); \
	done
	tar cvf - `cat .package` | bzip2 -9 > `cat .package`.tar.bz2
	@rm -rf `cat .package` .package

tar_agk:
	@echo $(PROJECT)_$(DATE) > .package
	@if [ -d "`cat .package`" ]; then \
	  echo ERROR: directory `cat .package` exits.  Stop. ;\
	else \
	  rm -rf `cat .package` ;\
	  mkdir `cat .package` `cat .package`/Makefiles ;\
	  for name in *.f90 *.fpp *.inc *.c makehead.awk test_os fortdep \
		AstroGK.in Makefile Makefile.depend Makefiles/* ;\
	  do ( ln $$name `cat .package`/$$name ; ) ;\
	  done ;\
	  echo tar cvf - `cat .package` \| bzip2 -9 \> `cat .package`.tar.bz2;\
	  tar cvf - `cat .package` | bzip2 -9 > `cat .package`.tar.bz2 ;\
	  rm -rf `cat .package` .package ;\
	fi

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
	@echo  USE_FFT is $(USE_FFT)
	@echo  USE_NETCDF is $(USE_NETCDF)
	@echo  USE_HDF5 is $(USE_HDF5)
	@echo  USE_C_INDEX is $(USE_C_INDEX)
	@echo  USE_NR_RAN is $(USE_NR_RAN)
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
	@echo PLATFORM_LINKS is $(PLATFORM_LINKS)

############################################################## PLATFORM LINKS
#
# If one of the platform-specific logical links doesn't exist, set them up:
# This is not used in agk any more
#
$(LINKS_$(PROJECT)):
	@$(MAKE) -C $(TOPDIR) --no-print-directory $(PLATFORM_LINKS)

#
# Platform-specific logical links:
#
cxt4: cxt4_$(PROJECT)
cxt4_gs2:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf ran_portable.f90 ran.f90
	ln -sf file_utils_portable.f90 file_utils.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf gs2_save_aix.f90 gs2_save.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_aix.f90 check.f90
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90 

c90: c90_$(PROJECT)
c90_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90 

t3e_shmem: t3e_shmem_$(PROJECT)
t3e_shmem_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

t3e: t3e_$(PROJECT)
t3e_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90 

t3e_fftw: t3e_fftw_$(PROJECT)
t3e_fftw_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

ibm: ibm_$(PROJECT)
ibm_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

origin: origin_$(PROJECT)
origin_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux: linux_$(PROJECT)
linux_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux_nomp: linux_nomp_$(PROJECT)
linux_nomp_gs2:
	ln -sf command_line_unix.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux_fuj: linux_fuj_$(PROJECT)
linux_fuj_gs2:
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
	ln -sf gs2_transforms_stub.f90 gs2_transforms.f90
	ln -sf fft_work_stub.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux_lf95: linux_lf95_$(PROJECT)
linux_lf95_gs2:
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
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux_abs: linux_abs_$(PROJECT)
linux_abs_gs2:
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
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

alpha: alpha_$(PROJECT)
alpha_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

alpha_nag: alpha_nag_$(PROJECT)
alpha_nag_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

linux_alpha: linux_alpha_$(PROJECT)
linux_alpha_gs2:
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

darwin: darwin_$(PROJECT)
darwin_gs2:
	ln -sf command_line_darwin.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -sf check_portable.f90 check.f90
	ln -sf ran_portable.f90 ran.f90
	ln -sf gs2_save_fast.f90 gs2_save.f90
	ln -sf file_utils_gnu.f90 file_utils.f90
	ln -sf gs2_transforms_fftw.f90 gs2_transforms.f90
	ln -sf fft_work_fftw.f90 fft_work.f90
	ln -sf gs2_dist_io_stub.f90 gs2_dist_io.f90
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

jacquard: jacquard_$(PROJECT)
jacquard_gs2:
	ln -sf command_line_unix.f90 command_line.f90
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
	ln -sf $(UTILS)/mds_io_stub.f90 mds.f90

unlink:	unlink_$(PROJECT)
unlink_gs2:
	-rm -f $(TOPDIR)/$(PROJECT)_dist_io.f90
	-rm -f $(TOPDIR)/$(PROJECT)_save.f90
	-rm -f $(TOPDIR)/$(PROJECT)_transforms.f90
	-rm -f $(TOPDIR)/command_line.f90
	-rm -f $(TOPDIR)/mp.f90
	-rm -f $(TOPDIR)/shmem.f90
	-rm -f $(TOPDIR)/prof.f90
	-rm -f $(TOPDIR)/redistribute.f90
	-rm -f $(TOPDIR)/ran.f90
	-rm -f $(TOPDIR)/fft_work.f90
	-rm -f $(TOPDIR)/check.f90 
	-rm -f $(TOPDIR)/file_utils.f90
	-rm -f $(TOPDIR)/mds.f90

unlink_agk:
	-rm -f $(TOPDIR)/$(PROJECT)_dist_io.f90
	-rm -f $(TOPDIR)/$(PROJECT)_io.f90
	-rm -f $(TOPDIR)/$(PROJECT)_layouts.f90
	-rm -f $(TOPDIR)/$(PROJECT)_layouts_type.h
	-rm -f $(TOPDIR)/$(PROJECT)_mem.f90
	-rm -f $(TOPDIR)/$(PROJECT)_save.f90
	-rm -f $(TOPDIR)/$(PROJECT)_transforms.f90
	-rm -f $(TOPDIR)/command_line.f90
	-rm -f $(TOPDIR)/constants.f90
	-rm -f $(TOPDIR)/fft_work.f90
	-rm -f $(TOPDIR)/file_utils.f90
	-rm -f $(TOPDIR)/mp.f90
	-rm -f $(TOPDIR)/ran.f90
