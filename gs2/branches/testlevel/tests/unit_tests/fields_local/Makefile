
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_fields_local

all: test_fields_local 
	$(call RUNCOMMAND,test_fields_local,test_fields_local.in)

#$(info "TEST_MODS=$(TEST_MODS)")

test_fields_local.o:  ${TEST_MODS} 

test_fields_local: test_fields_local.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

