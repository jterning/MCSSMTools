.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

aLibOBJ = fort.o hMasses.o nmssm.o

MDIR = $(subst /lib*?,,$(CURDIR)*?)

ifeq ($(SONAME),) 
  target = $(CURDIR)/mLib.$(SO)
else
  target = $(CURDIR)/mLib.$(SO) $(SONAME) $(CURDIR)/mLib.$(SO)
endif

all:  aLib.a mLib.$(SO)

aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f aLibOBJ

mLibSRC =  $(CALCHEP)/lib/model_aux.so

mLib.$(SO): $(mLibSRC)
override CFLAGS += $(CFLAGSD)
	$(CC) $(SHARED) $(CFLAGS) -o $(target) $(mLibSRC) $(MDIR)/../CalcHEP_src/lib/sqme_aux.$(SO) -lm 

clean: 
	rm -f *.o *.a *.so *.dll *.exe  so_locations 
