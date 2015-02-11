
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_theta_grid ogyropsi.dat EQDSK

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_theta_grid ogyropsi.dat EQDSK
	$(call RUNCOMMAND,test_theta_grid,)

endif

ogyropsi.dat EQDSK:
	tar -xjf data.tar.bz2

test_theta_grid.o:  ${TEST_MODS} 

test_theta_grid: test_theta_grid.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

