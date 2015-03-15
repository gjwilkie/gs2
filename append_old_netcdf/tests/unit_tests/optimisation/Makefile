
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_optimisation

all: test_optimisation 
	$(call RUNCOMMAND,test_optimisation,test_optimisation.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_optimisation.o:  ${TEST_MODS} 

test_optimisation: test_optimisation.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

