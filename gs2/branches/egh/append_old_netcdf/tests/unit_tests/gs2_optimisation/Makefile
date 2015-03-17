
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_optimisation

all: test_gs2_optimisation 
	$(call RUNCOMMAND,test_gs2_optimisation,test_gs2_optimisation.in)
	$(call RUNCOMMAND,../../../gs2,test_gs2_optimisation_run.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_optimisation.o:  ${TEST_MODS} 

test_gs2_optimisation: test_gs2_optimisation.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

