F90=f90

FFLAGS=-I$$MPTDIR/include -pgeo -putils 
#FFLAGS=-I$$MPTDIR/include -Iutils -64
#FFLAGS=-a taskcommon -I$$MPTDIR/include -pgeo -putils

FLIBS=-lmpi -lnetcdf 
#FLIBS=-L$$PAL_ROOT/lib64 -lVT -lmpi

LDFLAGS=
#LDFLAGS=-64

#COMMAND_LINE=command_line_unix
COMMAND_LINE=command_line_posix
#COMMAND_LINE=command_line_nag

CRAYF90HACK=:
#CRAYF90HACK=echo

MP=mp
#MP=mp_stub

PROF=prof_none
#PROF=prof_vampir

GS2MOD=gs2_diagnostics.o dist_fn.o collisions.o dist_fn_arrays.o \
	fields.o fields_arrays.o fields_implicit.o fields_test.o init_g.o \
	le_grids.o species.o run_parameters.o gs2_layouts.o $(MP).o \
	kt_grids.o kt_grids_single.o kt_grids_range.o kt_grids_specified.o \
	kt_grids_box.o \
	theta_grid.o theta_grid_file.o \
	theta_grid_salpha.o theta_grid_eik.o theta_grid_gridgen.o \
	theta_grid_params.o text_options.o file_utils.o $(COMMAND_LINE).o \
	$(PROF).o
GS2LIBMOD=utils/splines.o #geo/geo.a
#GS2LIB=gridgen/gridgen.a -Lutils -lsplines
GS2LIB=gridgen/gridgen.a utils/fitpack.o geo/radstub.o geo/vdimstub.o 

gs2: gs2.o $(GS2MOD)
	$(F90) $(LDFLAGS) -o gs2 gs2.o `$(CRAYF90HACK) $(GS2MOD)` $(GS2LIB) $(FLIBS)
#	$(F90) $(LDFLAGS) -o gs2 gs2.o `$(CRAYF90HACK) $(GS2MOD) $(GS2LIBMOD)` $(GS2LIB) $(FLIBS)

.SUFFIXES:
.SUFFIXES: .f90 .o .f
.f90.o:
	$(F90) -c $(FFLAGS) $<
.f.o:
	$(F90) -c $(FFLAGS) $<

kt_grids_range.o: file_utils.o
kt_grids_box.o: file_utils.o theta_grid.o
gs2_diagnostics.o: file_utils.o kt_grids.o run_parameters.o species.o $(MP).o
gs2_diagnostics.o: fields.o dist_fn.o constants.o $(PROF).o
dist_fn.o: $(MP).o species.o theta_grid.o kt_grids.o le_grids.o
dist_fn.o: run_parameters.o init_g.o
dist_fn.o: collisions.o gs2_layouts.o file_utils.o dist_fn_arrays.o constants.o
dist_fn.o: $(PROF).o
init_g.o: $(MP).o species.o theta_grid.o kt_grids.o le_grids.o dist_fn_arrays.o
init_g.o: gs2_layouts.o
run_parameters.o: $(MP).o file_utils.o
le_grids.o: $(MP).o species.o theta_grid.o kt_grids.o file_utils.o
le_grids.o: gs2_layouts.o constants.o
species.o: $(MP).o file_utils.o text_options.o
gs2.o: $(MP).o file_utils.o fields.o run_parameters.o gs2_diagnostics.o
file_utils.o: $(COMMAND_LINE).o
theta_grid_salpha.o: theta_grid_params.o file_utils.o text_options.o
theta_grid_salpha.o: constants.o theta_grid_gridgen.o
theta_grid_eik.o: theta_grid_params.o theta_grid_gridgen.o file_utils.o
fields.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o $(MP).o
fields.o: file_utils.o dist_fn_arrays.o constants.o $(PROF).o text_options.o
fields.o: fields_arrays.o fields_implicit.o fields_test.o
fields_implicit.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o $(MP).o
fields_implicit.o: fields_arrays.o theta_grid.o kt_grids.o dist_fn_arrays.o
fields_implicit.o: run_parameters.o gs2_layouts.o mp.o $(PROF).o file_utils.o
fields_test.o: theta_grid.o kt_grids.o run_parameters.o dist_fn.o $(MP).o
fields_test.o: fields_arrays.o theta_grid.o kt_grids.o dist_fn_arrays.o
fields_test.o: run_parameters.o gs2_layouts.o mp.o $(PROF).o file_utils.o
kt_grids.o: $(MP).o file_utils.o text_options.o kt_grids_single.o theta_grid.o
kt_grids.o: kt_grids_range.o kt_grids_specified.o kt_grids_box.o
kt_grids_single.o: file_utils.o
theta_grid_gridgen.o: file_utils.o constants.o
theta_grid_file.o: file_utils.o
theta_grid.o: $(MP).o file_utils.o text_options.o theta_grid_eik.o
theta_grid.o: theta_grid_salpha.o theta_grid_file.o
gs2_layouts.o: $(MP).o
theta_grid_params.o: file_utils.o
collisions.o: $(MP).o species.o theta_grid.o kt_grids.o le_grids.o
collisions.o: run_parameters.o file_utils.o text_options.o dist_fn_arrays.o
collisions.o: $(PROF).o
kt_grids_specified.o: file_utils.o
