
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_theta_grid ogyropsi.dat EQDSK

all: test_theta_grid ogyropsi.dat EQDSK
	$(call RUNCOMMAND,test_theta_grid,)

ogyropsi.dat EQDSK:
	tar -xjf data.tar.bz2

test_theta_grid.o:  ${TEST_MODS} 

test_theta_grid: test_theta_grid.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

