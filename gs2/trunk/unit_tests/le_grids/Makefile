
sinclude ../Makefile.tests_common
all: test_le_grids 
	mpirun -np $(NTESTPROCS) test_le_grids test_le_grids.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


test_le_grids.o:  ${TEST_MODS} 

test_le_grids: test_le_grids.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

