  
.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake
endif

all:micromegas.a 

override CFLAGS += $(LOCK) -D_SO_=\"$(SO)\"  -DMAKE=\"$(MAKE)\"

special = 
OBJ= $(notdir $(subst .c,.o,$(filter-out $(special), $(wildcard *.c )))) faux.o
micromegas.a:micromegas.a($(OBJ)) 
	$(RANLIB) micromegas.a 
	rm -f *.o

clean:
	rm -f *.a  *.o 
