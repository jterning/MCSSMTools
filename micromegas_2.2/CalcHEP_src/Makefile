
.PHONY: all  COMPILE clean

ifeq ($(MAKECMDGOALS),clean)
clean :
	if(test -r FlagsForSh) then rm -f  FlagsForSh; fi
	./setPath " "
	rm -f  FlagsForMake  so_locations CMessage
	cd lib; rm -f *.*
	cd bin; rm -f Int calc events2tab s_calchep show_distr sum_distr plot_view make_VandP  *.exe
	rm -f c_source/*/*.o c_source/*/so_location
endif

all:FlagsForMake COMPILE


FlagsForMake:
	./getFlags

COMPILE:FlagsForMake
	./setPath $(CURDIR)
	$(MAKE) -C c_source/getmem 
	$(MAKE) -C c_source/chep_crt
	$(MAKE) -C c_source/service2
	$(MAKE) -C c_source/polynom  
	$(MAKE) -C c_source/symb
	$(MAKE) -C c_source/plot 
	$(MAKE) -C c_source/num 
	$(MAKE) -C c_source/tab
	$(MAKE) -C c_source/model_aux
#	@if(test -r CMessage) then cat CMessage; fi
