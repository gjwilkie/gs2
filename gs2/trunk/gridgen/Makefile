INTERPOLATE=../utils
FFLAGS= -I$(INTERPOLATE) -p$(INTERPOLATE)
FC=f90
F90=f90
F90FLAGS= 
RANLIB=ranlib

GRIDGEN=gridgen4.o smooth.o
.SUFFIXES: .f90 .f
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(F90) $(F90FLAGS) -c $<

all: gridgen.a

gridgen.a: $(GRIDGEN)
	$(AR) rc $@ $(GRIDGEN)
	$(RANLIB) $@

# most systems
#command_line=command_line_unix
# nag
#command_line=command_line_nag
# cray mpp systems
command_line=command_line_posix

# cray f90 and similar
rungridgen=rungridgen.o
# digital f90 and similar
#rungridgen=rungridgen.o file_utils.o $(command_line).o

rungridgen.o: rungridgen.f90 file_utils.o
file_utils.o: file_utils.f90 $(command_line).o
$(command_line).o: $(command_line).f90

rungridgen: $(rungridgen) gridgen.a
	$(F90) -o $@ $(rungridgen) gridgen.a -L$(INTERPOLATE) -lsplines
