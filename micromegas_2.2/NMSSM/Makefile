.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake 
endif

all:
	$(MAKE) -C work
	$(MAKE) -C lib

clean: 
	$(MAKE) -C lib  clean
	$(MAKE) -C work clean 
	rm -f calchep/tmp/* calchep/results/*
	../sources/cleanexe
