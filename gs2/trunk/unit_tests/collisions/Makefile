
sinclude ../Makefile.tests_common
all: test_collisions 
	$(TESTCOMMAND) test_collisions test_collisions.in
	$(TESTCOMMAND) test_collisions test_collisions_use_le_layout.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


test_collisions.o:  ${TEST_MODS} 

test_collisions: test_collisions.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

