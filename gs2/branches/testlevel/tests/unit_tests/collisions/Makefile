
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_collisions

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_collisions 
	$(call RUNCOMMAND,test_collisions,test_collisions.in)
	$(call RUNCOMMAND,test_collisions,test_collisions_use_le_layout.in)

endif

test_collisions.o:  ${TEST_MODS} 

test_collisions: test_collisions.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

