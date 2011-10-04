.PHONY: clean init

ifeq ($(MAKECMDGOALS),init)
init:
	$(MAKE) -C micromegas_2.2
	$(MAKE) -C sources
else 

ifeq ($(MAKECMDGOALS),clean)
clean:
	$(MAKE) -C micromegas_2.2 clean
	$(MAKE) -C sources clean
	$(MAKE) -C main clean
else

micrO = micromegas_2.2

ifeq ($(wildcard $(micrO)/CalcHEP_src/FlagsForMake),)
$(error Use '[g]make init' for initialization)
endif

all:
	$(MAKE) -C main

endif

endif
