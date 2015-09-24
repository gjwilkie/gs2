
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_job_manage

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_job_manage 
	$(call RUNCOMMAND,test_job_manage,test_job_manage.in)

endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_job_manage.o:  ${TEST_MODS} 

test_job_manage: test_job_manage.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

