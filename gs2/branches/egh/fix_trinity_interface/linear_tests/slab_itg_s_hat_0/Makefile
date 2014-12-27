
sinclude ../Makefile.tests_common

CLEANCOMMAND+= slab_itg_s_hat_0_low_res

all: slab_itg_s_hat_0_low_res 
	$(call RUNCOMMAND,slab_itg_s_hat_0_low_res,slab_itg_s_hat_0_low_res.in)
	#mpirun -np $(NTESTPROCS) slab_itg_s_hat_0_high_res slab_itg_s_hat_0_high_res.in

slab_itg_s_hat_0_low_res.o:  ${TEST_MODS} 

slab_itg_s_hat_0_low_res: slab_itg_s_hat_0_low_res.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

