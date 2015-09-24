
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_diagnostics_new

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifdef USE_NEW_DIAG
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_gs2_diagnostics_new 
	$(call RUNCOMMAND,test_gs2_diagnostics_new,test_gs2_diagnostics_new.in)

endif

else

all:

endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_diagnostics_new.o:  ${TEST_MODS} 

test_gs2_diagnostics_new: test_gs2_diagnostics_new.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

