
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_theta_grid

all: test_theta_grid 
	$(call RUNCOMMAND,test_theta_grid,test_theta_grid.in)

test_theta_grid.o:  ${TEST_MODS} 

test_theta_grid: test_theta_grid.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

