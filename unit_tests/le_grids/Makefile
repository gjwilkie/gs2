
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_le_grids

all: test_le_grids 
	$(call RUNCOMMAND,test_le_grids,test_le_grids.in)

test_le_grids.o:  ${TEST_MODS} 

test_le_grids: test_le_grids.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

