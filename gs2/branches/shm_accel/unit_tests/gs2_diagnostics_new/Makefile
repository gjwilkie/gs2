
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_diagnostics_new

all: test_gs2_diagnostics_new 
	$(call RUNCOMMAND,test_gs2_diagnostics_new,test_gs2_diagnostics_new.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_diagnostics_new.o:  ${TEST_MODS} 

test_gs2_diagnostics_new: test_gs2_diagnostics_new.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

