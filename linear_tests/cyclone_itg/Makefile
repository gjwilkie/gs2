
sinclude ../Makefile.tests_common

CLEANCOMMAND+= cyclone_itg_low_res

all: cyclone_itg_low_res 
	$(call RUNCOMMAND,cyclone_itg_low_res,cyclone_itg_low_res.in) implicit
	$(call RUNCOMMAND,cyclone_itg_low_res,cyclone_itg_low_res_fields_local.in) local
	#mpirun -np $(NTESTPROCS) cyclone_itg_high_res cyclone_itg_high_res.in

cyclone_itg_low_res.o:  ${TEST_MODS} 

cyclone_itg_low_res: cyclone_itg_low_res.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

