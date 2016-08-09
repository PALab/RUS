sources = rus_alloc.c \
	  rus_pars.c \
	  xindex.c \
	  rus_forward.c \
	  rus_inverse.c

objects = $(sources:.c=.o)

makes = $(sources:.c=.d)

default: $(objects)

-include $(makes)

CC=gcc

%.d: %.c
	@set -e; rm -f $@; \
         $(CC) -M $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
         rm -f $@.$$$$

clean:
	-rm $(objects) $(makes)
