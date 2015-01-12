
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_job_manage

all: test_job_manage 
	$(call RUNCOMMAND,test_job_manage,test_job_manage.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_job_manage.o:  ${TEST_MODS} 

test_job_manage: test_job_manage.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

