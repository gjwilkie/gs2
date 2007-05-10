#
# NOTE: The file test_os should be executable.
#
#################################################################### OVERVIEW
#  Makefile for the GS2 Gyrokinetic Stability code 
#  (requires GNU's gmake, SEE Makefile.basic for the basic version).  
#
#  Written by D. Ernst <dernst@pppl.gov> 3/99, with input from 
#  Greg Hammett and Bill Dorland 
#
#  LAST UPDATE: 4/8/99
################################################################# TEST RESULTS
# TEST RESULTS as of 4/99:
#
#  Cray T3E: runs, no problems
# 
#  SGI Origin 2000: runs, no problems

#  Linux, NAG f95: compiles, does not run
#   Makefile works.
#   gs2 crashes with I/O problem (Does grid.out exist in test directory? yes)
#   ball runs
#   
#     ./gs2 s1.in
#     Running on  1  processor
#     Run name: s1
#     could not open grid.out file (1)
#     Unit 51 is not connected
#     Program terminated by fatal I/O error
#     Aborted (core dumped)
# 
#  DEC Alpha: compiles, does not run  
#
#  (but Greenwald has working code for DEC ALPHA under VMS)
#
#     Makefile works. 
#     gs2 crashes with seg fault.
#     > "limit stacksize unlimited" does not help
#     Doesn't crash in totalview.
# 
#     Running on            1  processor
#     Run name: gs  # gs.in is a symlink to s1.in
#     forrtl: severe (139): array index out of bounds for index 1 (SIGTRAP)
#
#  DEC Alpha, f95:
#
#  Cray J90:
############################################################# WALL CLOCK TIMES
# WALL CLOCK TIMES for s1.in (7 modes)
#
#proc   T3E     Origin  Linux-266MMX  Alpha-500MHz  Alpha-f95-500MHz J90
# 1     34m     33:44m     
# 4     9:39m   11:07m
# 8     6m      6:00m
# 16    3:47m   4:12m
# 32    2:48m   
# 64    1:58m
# 128   
#
######################################################################## TODO
#  TODO:
#      -copy over netcdf_dummy.f90
#      -add to links where needed
#
############################################################ MAKE INSTRUCTIONS
#
#  First do "module load netcdf " (on Crays) or otherwise set up your
#  compilation/library environment.
#
#  Then type:
#
#    "make" to build an optimized version of GS2
#    "make debug=on" to build debug version
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
# NEW FEATURES:
#
#   "make test_make" - print values of variables
#   "make tar" - prepare a dated tar.gz file of src directory 
#   "make distclean" - does "make clean" + removes platform links
# 
####################################################################### CPU
#
# environment variable CPU is used to set platform-dependent make variables.
#
# If CPU is undefined; use the test_os script to figure it out:
ifeq ($(CPU),)
   CPU_tmp := $(shell ./test_os)
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
  ifneq ($(debug),on)
    FFLAGS =  -N80 -I/usr/local/include -M1110 -O vector3 -O inline4
    F90FLAGS =  -I/usr/local/include -M1110 -O vector3 -O inline4
    FLIBS = -lmpi -L/usr/local/lib  $$NETCDF 
  else
    FFLAGS =  -N80 -I/usr/local/include -M1110 -g -R abcs -e i
    F90FLAGS =     -I/usr/local/include -M1110 -g -R abcs -e i
    FLIBS   = -lmpi -L/usr/local/lib $$NETCDF -Wl"-D preset=inf" # -lmalloc 
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
  ifneq ($(debug),on)
    FFLAGS =  -N80 -I/usr/local/include -M1110 -O vector3 
    F90FLAGS =  -I/usr/local/include -M1110 -O vector3  
    FLIBS = -L/usr/local/lib $$NETCDF
  else
    FFLAGS =  -N80 -I/usr/local/include -M1110 -g -R abcs -e i
    F90FLAGS =     -I/usr/local/include -M1110 -g -R abcs -e i
    FLIBS   = -L/usr/local/lib  $$NETCDF  -Wl"-D preset=inf" # -lmalloc
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
  ifneq ($(debug),on)
    FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 -O -TARG:platform=ip27
    F90FLAGS =  -I/usr/pppl/include -64 -r8 -O -TARG:platform=ip27
# Other options tried (no more than 15% speedup):
#    FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
#    F90FLAGS =  -I/usr/pppl/include -64 -r8 -Ofast=ip27 \
#	-TARG:platform=ip27 -OPT:IEEE_arithmetic=3 -lfastm
  else
    FFLAGS =  -col80 -I/usr/pppl/include -64 -r8 -g \
	-DEBUG:div_check=3:trap_uninitialized=on
    F90FLAGS =  -I/usr/pppl/include -64 -r8 -g \
	-DEBUG:div_check=3:trap_uninitialized=on
#       -DEBUG:div_check=3:subscript_check=on:trap_uninitialized=on
# should be the full debug options, but then the compiler bombs on some
# of the routines (even the new beta version 7.3 of the compiler).
  endif
endif

# DEC alpha options:
ifeq ($(CPU),ALPHA)
  FC = f95
  FLIBS  = -L/usr/local/lib -lnagdx -lnetcdf -ldxml
  NAG_DP = nag_dp.o
  PLATFORM_LINKS = alpha
  ifeq ($(debug),on)
    FFLAGS =  -extend_source -I/usr/local/include -r8 -g \
	-assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source -align dcommons
    F90FLAGS =               -I/usr/local/include -r8 -g \
	-assume dummy_aliases -check bounds -check overflow \
	-warn argument_checking -warn truncated_source -align dcommons
  else
     FFLAGS =  -extend_source -I/usr/local/include -r8 -O -arch host \
	-align dcommons -math_library fast
     F90FLAGS =               -I/usr/local/include -r8 -O -arch host \
	-align dcommons -math_library fast 
  endif
endif

# options for NAG f95 on a DEC alpha
ifeq ($(CPU),ALPHA_NAG)
  FC = /usr/local/bin/f95
  FLIBS = -L/usr/local/lib -lnetcdf
  PLATFORM_LINKS = alpha_nag
  ifeq ($(debug),on)
    FFLAGS = -C -132 -I/usr/local/include -r8 -g90 -dusty
    F90FLAGS = -C -hpf  -I/usr/local/include -r8 -g90 -dusty
  else
    FFLAGS =  -C -132 -I/usr/local/include -r8 -O -dusty
    F90FLAGS =  -C -hpf -I/usr/local/include -r8 -O -dusty
  endif
endif

# options for Linux with NAG f95:
ifeq ($(CPU),LINUX)
  FC = f95
  FLIBS = -L/usr/local/lib -L/usr/lib -lnetcdf
  PLATFORM_LINKS = linux
  ifeq ($(debug),on)
    FFLAGS =  -C -132 -I/usr/local/include -r8 -g90 -gline -P -dusty
    F90FLAGS = -C  -hpf -I/usr/local/include -r8 -g90 -gline -P -dusty
  else
    FFLAGS =  -C -132 -I/usr/local/include -r8 -O -dusty
    F90FLAGS = -C -hpf -I/usr/local/include -r8 -O -dusty
  endif
endif

# 
# SUN or DEC alpha switches?:
# -eA -g -dalign -I/usr/local/include 
# -ef #-Rabc -m0

########################################################## MODULE DECLARATIONS

# Cray compiler automatically searches for module *.o files in present
# directory, but NAG compiler doesn't, so have to explicitly list them:

stubs = vdimstub.o radstub.o
eqmod = veq.o geq.o eeq.o leq.o peq.o 

GS2MOD=gs2_diagnostics.o dist_fn.o collisions.o dist_fn_arrays.o \
	fields.o fields_arrays.o fields_implicit.o fields_test.o init_g.o \
	le_grids.o species.o run_parameters.o geometry.o \
        gs2_layouts.o \
	kt_grids.o kt_grids_single.o kt_grids_range.o kt_grids_specified.o \
	kt_grids_box.o kt_grids_xbox.o $(eqmod) \
	theta_grid.o theta_grid_file.o gs2_reinit.o \
	theta_grid_salpha.o theta_grid_eik.o theta_grid_gridgen.o \
	theta_grid_params.o gs2_save.o gs2_time.o redistribute.o \
        text_options.o file_utils.o ran.o \
	command_line.o prof.o splines.o mp.o shmem.o

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

all: eiktest ball rungridgen gs2

eiktest: et.o fitpack.o  abort.o
	$(FC) $(FFLAGS) -o  eiktest et.o fitpack.o geometry.o $(eqmod) \
        splines.o abort.o $(stubs) $(FLIBS)

ball: ball.o fitpack.o abort.o
	$(FC) $(FFLAGS) -o ball ball.o  fitpack.o geometry.o $(eqmod) \
        splines.o abort.o $(stubs) $(FLIBS)

rungridgen: rungridgen.o gridgen4.o smooth.o fitpack.o \
	file_utils.o command_line.o abort.o
	$(FC) $(FFLAGS) -o rungridgen rungridgen.o \
	gridgen4.o smooth.o fitpack.o file_utils.o \
	command_line.o abort.o $(FLIBS)

gs2: gs2.o $(GS2MOD) fitpack.o gridgen4.o abort.o
	$(FC) $(F90FLAGS) -o gs2 gs2.o gridgen4.o fitpack.o abort.o \
	$(stubs) $(GS2MOD) $(FLIBS)


################################################################# DEPENDENCIES

kt_grids_range.o: file_utils.o
kt_grids_box.o: file_utils.o theta_grid.o
kt_grids_xbox.o: file_utils.o theta_grid.o
gs2_diagnostics.o: file_utils.o kt_grids.o run_parameters.o species.o mp.o 
gs2_diagnostics.o: fields.o dist_fn.o constants.o prof.o gs2_save.o gs2_time.o
dist_fn.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o
dist_fn.o: run_parameters.o init_g.o
dist_fn.o: gs2_layouts.o file_utils.o dist_fn_arrays.o constants.o
dist_fn.o: collisions.o 
dist_fn.o: prof.o gs2_time.o
init_g.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o dist_fn_arrays.o
init_g.o: gs2_layouts.o gs2_save.o fields_arrays.o ran.o
run_parameters.o: mp.o file_utils.o
le_grids.o: mp.o species.o theta_grid.o kt_grids.o file_utils.o redistribute.o
le_grids.o: gs2_layouts.o constants.o 
species.o: mp.o file_utils.o text_options.o
gs2.o: mp.o file_utils.o fields.o run_parameters.o gs2_diagnostics.o 
gs2.o: gs2_reinit.o gs2_time.o init_g.o
file_utils.o: command_line.o
theta_grid_salpha.o: theta_grid_params.o file_utils.o text_options.o
theta_grid_salpha.o: constants.o theta_grid_gridgen.o
theta_grid_eik.o: theta_grid_params.o theta_grid_gridgen.o file_utils.o
fields.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields.o: file_utils.o dist_fn_arrays.o constants.o prof.o text_options.o
fields.o: fields_arrays.o fields_implicit.o fields_test.o 
fields_implicit.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields_implicit.o: fields_arrays.o theta_grid.o kt_grids.o dist_fn_arrays.o
fields_implicit.o: run_parameters.o gs2_layouts.o mp.o prof.o file_utils.o
fields_implicit.o: gs2_save.o
fields_test.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o mp.o
fields_test.o: fields_arrays.o theta_grid.o kt_grids.o dist_fn_arrays.o
fields_test.o: run_parameters.o gs2_layouts.o mp.o prof.o file_utils.o
kt_grids.o: mp.o file_utils.o text_options.o kt_grids_single.o theta_grid.o
kt_grids.o: kt_grids_range.o kt_grids_specified.o kt_grids_box.o
kt_grids.o: kt_grids_xbox.o
kt_grids_single.o: file_utils.o
theta_grid_gridgen.o: file_utils.o constants.o splines.o
theta_grid_file.o: file_utils.o
theta_grid.o: mp.o file_utils.o text_options.o theta_grid_eik.o
theta_grid.o: theta_grid_salpha.o theta_grid_file.o
gs2_layouts.o: mp.o
theta_grid_params.o: file_utils.o
collisions.o: mp.o species.o theta_grid.o kt_grids.o le_grids.o
collisions.o: run_parameters.o file_utils.o text_options.o dist_fn_arrays.o
collisions.o: prof.o shmem.o redistribute.o
kt_grids_specified.o: file_utils.o
theta_grid_eik.o: theta_grid_params.o theta_grid_gridgen.o file_utils.o geometry.o
gs2_save.o: theta_grid.o gs2_layouts.o mp.o fields_arrays.o kt_grids.o dist_fn_arrays.o \
	gs2_time.o
gs2_reinit.o: collisions.o mp.o \
	dist_fn.o fields.o fields_implicit.o fields_test.o init_g.o run_parameters.o \
	gs2_save.o dist_fn_arrays.o fields_arrays.o 
redistribute.o: mp.o


et.o: geometry.o 
ball.o: geometry.o
geometry.o: veq.o geq.o eeq.o leq.o peq.o radstub.o
veq.o: splines.o vdimstub.o
geq.o: splines.o 
eeq.o: splines.o

rungridgen.o: file_utils.o

############################################################## MORE DIRECTIVES
clean:
	rm -f *.o *~ *.mod *.g90 core */core

distclean:	unlink clean

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

%.o : %.f
	$(FC) $(FFLAGS) -c $<

test_make:
	@echo FC is $(FC)
	@echo FFLAGS is $(FFLAGS)
	@echo F90FLAGS is $(F90FLAGS)
	@echo debug is $(debug)
	@echo CPU is $(CPU)

# If one of the platform-specific logical links doesn't exist, set them up:

command_line.f90 mp.f90 shmem.f90 prof.f90:
	$(MAKE) $(PLATFORM_LINKS)

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
	ln -s ran_cray.f90 ran.f90

t3e_shmem:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi.f90 mp.f90
	ln -sf shmem_cray.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_shmem.f90 redistribute.f90
	ln -s ran_cray.f90 ran.f90

t3e:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -s ran_cray.f90 ran.f90

origin:
	ln -sf command_line_posix.f90 command_line.f90
	ln -sf mp_mpi_r8.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -s ran_cray.f90 ran.f90

linux:
	ln -sf command_line_nag.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -s ran_portable.f90 ran.f90

alpha:
	ln -sf command_line_alpha.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -s ran_portable.f90 ran.f90

alpha_nag:
	ln -sf command_line_nag.f90 command_line.f90
	ln -sf mp_stub.f90 mp.f90
	ln -sf shmem_stub.f90 shmem.f90
	ln -sf prof_none.f90 prof.f90
	ln -sf redistribute_mpi.f90 redistribute.f90
	ln -s ran_portable.f90 ran.f90

unlink:
	rm -f command_line.f90
	rm -f mp.f90
	rm -f shmem.f90
	rm -f prof.f90
	rm -f redistribute.f90
	rm -f ran.f90
