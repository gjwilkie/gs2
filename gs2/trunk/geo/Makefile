# Generated automatically from Makefile.in by configure.
#
# Makefile for geometry.a
#
FC=f90
FLIBS= -lnetcdf
FFLAGS=  -N80  -I/usr/local/include -I.. -I.. -p../utils 
RANLIB=ranlib
NTDIR=..
UTILDIR=../utils

.SUFFIXES: .f90 .f

.f90.o:
	f90 $(FFLAGS) -c $<

.f.o:
	f90 $(FFLAGS) -c $<

stubs = radstub.o vdimstub.o

geo_obj = geometry.o geq.o veq.o eeq.o peq.o leq.o 

geo.a: $(geo_obj) $(stubs)
	$(AR) rc $@ $(geo_obj) $(stubs)
	$(RANLIB) $@

all:	geo.a

clean:
	rm -f *.o *~ geometry.a  *.tmp *.trace core eiktest

distclean: clean
	$(RM) nt.par 
	$(RM) config.cache config.log config.status Makefile

configure: configure.in
	autoconf --local=config

geometry.o : veq.o geq.o eeq.o leq.o peq.o 

ts : testsmoothf.o smoothf.o
	$(FC) $(FFLAGS) -o ts testsmoothf.o smoothf.o 

ball: ball.o geo.a ../utils/fitpack.o $(stubs)
	$(FC) $(FFLAGS) -o ball ball.o $(stubs) \
	../utils/fitpack.o $(FLIBS) 

eiktest: mods et.o geo.a ../utils/fitpack.o $(stubs)
	$(FC) -o eiktest et.o ../utils/fitpack.o $(stubs) $(FLIBS)

et.o: geometry.o 

ball.o: geometry.o 

peq.o geq.o: $(UTILDIR)/splines.a $(UTILDIR)/netcdf.a 
veq.o: $(UTILDIR)/splines.a
eeq.o: $(UTILDIR)/splines.a

MODULES= $(UTILDIR)

mods: 
	for dir in $(MODULES); do (cd $$dir; $(MAKE)); done



