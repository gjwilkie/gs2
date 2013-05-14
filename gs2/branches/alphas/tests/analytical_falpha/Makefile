
.SUFFIXES: .fpp .f90 .c .o

F90FLAGS += -I ../../
LDFLAGS += -I ../../
CLAGS += -I ../../

.f90.o: 
	$(FC) $(F90FLAGS) -c $<
.fpp.f90:
	$(CPP) $(CPPFLAGS) $< $@
.c.o:
	$(CC) $(CFLAGS) -c $<

all: test_analytical_falpha

test_analytical_falpha.o: ../../analytical_falpha.o 

test_analytical_falpha: test_analytical_falpha.o ../../analytical_falpha.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
	./test_analytical_falpha

