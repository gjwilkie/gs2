
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_reinit

all: test_gs2_reinit 
	$(call RUNCOMMAND,test_gs2_reinit,test_gs2_reinit.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_reinit.o:  ${TEST_MODS} 

test_gs2_reinit: test_gs2_reinit.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

