.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../micromegas_2.2/CalcHEP_src/FlagsForMake
endif

micrO = ../micromegas_2.2

LIB = ../sources/nmhdecay.a \
      $(micrO)/sources/micromegas.a  $(CALCHEP)/lib/sqme_aux.so  \
      $(micrO)/NMSSM/lib/aLib.a   $(CALCHEP)/lib/model_aux.so  \
      $(micrO)/NMSSM/work/work_aux.a    


EXEC = nmhdecay nmspec nmgmsb nmhdecay_grid nmhdecay_rand nmspec_grid nmspec_rand nmgmsb_grid nmgmsb_rand

all: $(EXEC)

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

nmhdecay:	nmhdecay.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmhdecay_grid:	nmhdecay_grid.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmhdecay_rand:	nmhdecay_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmspec: 	nmspec.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmspec_grid: 	nmspec_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmspec_rand: 	nmspec_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmgmsb: 	nmgmsb.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmgmsb_grid: 	nmgmsb_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

nmgmsb_rand: 	nmgmsb_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL)

clean:
	rm -f *.o $(EXEC)
