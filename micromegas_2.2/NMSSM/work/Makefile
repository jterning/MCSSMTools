
.PHONY:  clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

all: work_aux.a newProcess

workOBJ =  VandP.o  path.o
override CFLAGS += -I../../CalcHEP_src/include

work_aux.a: work_aux.a($(workOBJ))
	$(RANLIB) work_aux.a 
	rm -f $(workOBJ)

VandP.c: models/vars1.mdl models/func1.mdl models/prtcls1.mdl
	../../CalcHEP_src/bin/make_VandP models 1 

path.c:
	echo char \* WORK=\"`pwd`\"\; >  path.c

newProcess: ../../sources/newProcess
	cp  ../../sources/newProcess .

clean:
	rm -f *.a *.o  lock_ LOCK  results/* tmp/* models/*2.mdl so_generated/* lock_
	rm -f  path.c VandP.c newProcess
