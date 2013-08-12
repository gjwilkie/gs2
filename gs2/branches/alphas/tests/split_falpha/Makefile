
sinclude ../Makefile.tests_common
all: test_split_falpha
	mpirun -np $(NTESTPROCS) $^
	


test_split_falpha.o: ../../split_falpha.o $(TEST_MODS)

test_split_falpha: test_split_falpha.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

