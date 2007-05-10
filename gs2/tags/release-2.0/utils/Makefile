#
# Makefile for utils.a
#
FC=f90
#FFLAGS=@FFLAGS@
FFLAGS = 
RANLIB=ranlib
#NTDIR=@NTDIR@

.SUFFIXES: .f90 .f

.f90.o:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

all:	utils.a splines.a netcdf.a 

splines_obj = splines.o fitpack.o 

splines.a: $(splines_obj) 
	$(AR) rc $@ $(splines_obj)
	$(RANLIB) $@

netcdf.a: netcdf.o 
	$(AR) rc $@ netcdf.o
	$(RANLIB) $@

utils_obj = aborter.o lnblank.o gheat.o #@LAPACK@ @NAG@

utils.a: $(utils_obj)
	$(AR) rc $@ $(utils_obj)
	$(RANLIB) $@

tix2nt_obj = tix2nt.o u2read.o

tix2nt: $(tix2nt_obj)
	$(FC) $(FFLAGS) -o tix2nt $(tix2nt_obj) utils.a $(FLIBS)

clean:
	rm -f *.o *~ *.a  *.tmp *.trace core

distclean: clean
	$(RM) config.cache config.log config.status Makefile
	$(RM) netcdf.inc

configure: configure.in
	autoconf --local=config
