
sinclude ../Makefile.tests_common
all: test_general_f0
	./test_general_f0

MODS=../../unit_tests.o ../../general_f0.o ../../mp.o ../../spl.o ../../file_utils.o 

MODS=${echo gs2_mod | sed s/\(.+\)/..\/..\/\1/}
MODS=${gs2_mod}

test_general_f0.o:  ${MODS} 

test_general_f0: test_general_f0.o  ${MODS}
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

