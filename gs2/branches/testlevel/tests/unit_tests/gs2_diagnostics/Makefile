
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_diagnostics

all: test_gs2_diagnostics 
	$(call RUNCOMMAND,test_gs2_diagnostics,test_gs2_diagnostics.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_diagnostics.o:  ${TEST_MODS} 

test_gs2_diagnostics: test_gs2_diagnostics.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

