
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_collisions

all: test_collisions 
	$(call RUNCOMMAND,test_collisions,test_collisions.in)
	$(call RUNCOMMAND,test_collisions,test_collisions_use_le_layout.in)

test_collisions.o:  ${TEST_MODS} 

test_collisions: test_collisions.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

