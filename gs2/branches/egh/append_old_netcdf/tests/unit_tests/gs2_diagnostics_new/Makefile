
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_diagnostics_new old_diagnostics/.test_* old_diagnostics/*.{o,error,out*,vres,*mom*,eigenfunc,fields,g,heat*,jext,lpc,parity,phase} averages.dat tmpdata.dat

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else

# Here we are essentially running two tests. First we run the old
# diagnostics. Then we run the new diagnostics and check the results
# are the same. The we run the new diagnostics for half the number of
# steps. Then we restart and run for another half, copying the file
# from the first half and appending to it, and check that the result
# is still unchanged.
all: test_gs2_diagnostics_new 
	$(call RUNCOMMAND,test_gs2_diagnostics_new,old_diagnostics/test_gs2_diagnostics_new.in)
	$(call RUNCOMMAND,test_gs2_diagnostics_new,test_gs2_diagnostics_new.in)
	rm -f test_gs2_diagnostics_new_start.out.nc
	$(call RUNCOMMAND,test_gs2_diagnostics_new,test_gs2_diagnostics_new_start.in)
	cp test_gs2_diagnostics_new_start.out.nc test_gs2_diagnostics_new_append.out.nc; echo
	$(call RUNCOMMAND,test_gs2_diagnostics_new,test_gs2_diagnostics_new_append.in)

endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_diagnostics_new.o:  ${TEST_MODS} 

test_gs2_diagnostics_new: test_gs2_diagnostics_new.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

