.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../micromegas_2.2/CalcHEP_src/FlagsForMake
endif

all:	nmhdecay.a 

SRC = $(wildcard *.f)
OBJ = $(SRC:.f=.o)

nmhdecay.a:nmhdecay.a($(OBJ)) 
	$(RANLIB) nmhdecay.a 
	rm -f *.o

clean:
	rm -f *.a  *.o
