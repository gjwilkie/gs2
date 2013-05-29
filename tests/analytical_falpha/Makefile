
sinclude ../Makefile.tests_common
all: test_analytical_falpha
	mpirun -np $(NTESTPROCS) $^
	


test_analytical_falpha.o: $(TEST_MODS)

test_analytical_falpha: test_analytical_falpha.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

