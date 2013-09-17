
sinclude ../Makefile.tests_common
all: wstar_units 
	$(TESTCOMMAND) wstar_units wstar_units.list
	#mpirun -np $(NTESTPROCS) cyclone_itg_high_res cyclone_itg_high_res.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


wstar_units.o:  ${TEST_MODS} 

wstar_units: wstar_units.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

