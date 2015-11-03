sinclude ../Makefile.tests_common
all: time_ffts 
	$(TESTCOMMAND) time_ffts time_ffts.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o



time_ffts.o:  ${TEST_MODS} 

time_ffts: time_ffts.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

