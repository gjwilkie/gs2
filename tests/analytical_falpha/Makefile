
sinclude ../Makefile.tests_common
all: test_analytical_falpha
	./test_analytical_falpha


test_analytical_falpha.o: ../../analytical_falpha.o ../../unit_tests.o

test_analytical_falpha: test_analytical_falpha.o ../../analytical_falpha.o ../../unit_tests.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

