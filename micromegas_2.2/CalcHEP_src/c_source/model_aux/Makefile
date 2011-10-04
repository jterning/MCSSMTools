.PHONY: all

include ../../FlagsForMake

all: $(CALCHEP)/lib/model_aux.$(SO)


ifneq ($(SONAME),)
ifeq ($(CC),gcc)
   soname_model_aux :=  -Wl,$(SONAME),$(CALCHEP)/lib/model_aux.$(SO)
else
   soname_model_aux :=   $(SONAME)   $(CALCHEP)/lib/model_aux.$(SO)
endif
else
   soname_model_aux :=  
endif

$(CALCHEP)/lib/model_aux.$(SO): alpha_s.c  rdiag.c  SLHAreader.c
	$(CC) $(CFLAGS) -Iinclude $(SHARED) -o $@ $(soname_model_aux)  $^  $(CALCHEP)/lib/sqme_aux.$(SO) -lm
 
