
sinclude ../Makefile.tests_common

CLEANCOMMAND+= wstar_units

all: wstar_units 
	$(call RUNCOMMAND,wstar_units,wstar_units.list)

wstar_units.o:  ${TEST_MODS} 

wstar_units: wstar_units.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

