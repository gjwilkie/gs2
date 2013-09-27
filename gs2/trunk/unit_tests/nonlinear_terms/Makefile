
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_nonlinear_terms

all: test_nonlinear_terms 
	$(call RUNCOMMAND,test_nonlinear_terms,test_nonlinear_terms.in)

test_nonlinear_terms.o:  ${TEST_MODS} 

test_nonlinear_terms: test_nonlinear_terms.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

