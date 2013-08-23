
# This command adds '../../' to all the gs2 modules and gets rid of gs2.o
TEST_MODS=$(shell echo '$(gs2_mod)' | sed -e "s/gs2.o//" |  sed -e "s/\(\S\+\)/..\/..\/\1/g" )
TEST_MODS+= ../../unit_tests.o ../../functional_tests.o

# All subfolders unless otherwise specified
TESTS ?= $(wildcard */.) 

     
.PHONY: all $(TESTS)

all: $(TESTS)

export
$(TESTS): 
	$(MAKE) -C $@
