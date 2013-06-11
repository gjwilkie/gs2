
sinclude ../Makefile.tests_common
all: test_general_f0
	mpirun -np $(NTESTPROCS) test_general_f0 test_general_f0.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o ../../run_parameters.o


test_general_f0.o:  ${TEST_MODS} 

test_general_f0: test_general_f0.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

