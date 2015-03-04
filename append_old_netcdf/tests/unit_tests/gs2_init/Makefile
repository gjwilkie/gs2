
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_init

all: test_gs2_init 
	$(call RUNCOMMAND,test_gs2_init,test_gs2_init.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_init.o:  ${TEST_MODS} 

test_gs2_init: test_gs2_init.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

