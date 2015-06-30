
sinclude ../Makefile.tests_common

CLEANCOMMAND+= cyclone_itg_low_res

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: cyclone_itg_low_res cyclone_itg_replay
	#$(call RUNCOMMAND,cyclone_itg_low_res,cyclone_itg_low_res.in,implicit)
	#$(call RUNCOMMAND,cyclone_itg_low_res,cyclone_itg_low_res_fields_local.in,local)
	cp cyclone_itg_low_res_fields_local.out.nc cyclone_itg_replay.out.nc
	$(call RUNCOMMAND,cyclone_itg_replay,cyclone_itg_replay.in,)
endif

cyclone_itg_low_res.o cyclone_itg_replay.o:  ${TEST_MODS} 

cyclone_itg_low_res: cyclone_itg_low_res.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

cyclone_itg_replay: cyclone_itg_replay.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

