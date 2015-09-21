
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_le_grids

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_le_grids 
	$(call RUNCOMMAND,test_le_grids,test_le_grids.in)

endif

test_le_grids.o:  ${TEST_MODS} 

test_le_grids: test_le_grids.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

