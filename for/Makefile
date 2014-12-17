# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config

B = $(HOME)/bin

D = $L/libcwp.a $L/libpar.a $L/libsu.a

OPTC = -O0 -Wall -ansi -pedantic  -DHAVE_UNISTD_H

#LFLAGS= -lm  -L$L -lsu -lpar -lcwp -lblas -llapack -lg2c -lI77 -lF77
LFLAGS= -lm  -L$L -lsu -lpar -lcwp -lblas -llapack #-lg2c

PROGS =			\
	$B/formod_aniso




INSTALL	:	$(PROGS)
	touch $@


$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@
	@chmod 755 $@
	@echo $(@F) installed in $B

remake	:
	-touch *.c 
	-rm -f $(PROGS)
	$(MAKE) 
clean:
	rm -f a.out junk* JUNK* core
