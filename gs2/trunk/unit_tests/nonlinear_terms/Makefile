
sinclude ../Makefile.tests_common
all: test_nonlinear_terms 
	$(TESTCOMMAND) test_nonlinear_terms test_nonlinear_terms.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


test_nonlinear_terms.o:  ${TEST_MODS} 

test_nonlinear_terms: test_nonlinear_terms.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

