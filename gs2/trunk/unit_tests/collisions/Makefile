
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_collisions

all: test_collisions 
	$(TESTCOMMAND) test_collisions test_collisions.in
	$(TESTCOMMAND) test_collisions test_collisions_use_le_layout.in

test_collisions.o:  ${TEST_MODS} 

test_collisions: test_collisions.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

