
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_nonlinear_terms

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_nonlinear_terms 
	$(call RUNCOMMAND,test_nonlinear_terms,test_nonlinear_terms.in)

endif

test_nonlinear_terms.o:  ${TEST_MODS} 

test_nonlinear_terms: test_nonlinear_terms.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

