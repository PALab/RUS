
# Makefile for ...su/main

include $(CWPROOT)/src/Makefile.config

B = $(HOME)/bin

D = $L/libcwp.a $L/libpar.a $L/libsu.a

OPTC =  -g -Wall -ansi -pedantic  -DHAVE_UNISTD_H

LFLAGS= -L$L -lpar -lcwp -lsu -lm -lblas -llapack #-lg2c 
# -lpar -lsu -lcwp   		Get from CWP su download 
# -lm 				Standard C math package
# -lblas -llapack  -lg2c	Get from RPMfind.com, linear alg.

PROGS =			\
	$B/RUS-inverse

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
