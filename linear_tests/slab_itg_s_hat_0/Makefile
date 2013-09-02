
sinclude ../Makefile.tests_common
all: slab_itg_s_hat_0_low_res 
	$(TESTCOMMAND) slab_itg_s_hat_0_low_res slab_itg_s_hat_0_low_res.in
	#mpirun -np $(NTESTPROCS) slab_itg_s_hat_0_high_res slab_itg_s_hat_0_high_res.in



# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


slab_itg_s_hat_0_low_res.o:  ${TEST_MODS} 

slab_itg_s_hat_0_low_res: slab_itg_s_hat_0_low_res.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

