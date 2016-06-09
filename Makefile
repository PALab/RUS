CC = gcc

LFLAGS = -lm -llapack

.PHONY : default clean install uninstall \
    forward_example inverse_example objects

default: rus_forward rus_inverse

objects:
	$(MAKE) -C src

install: rus_forward rus_inverse
	cp rus_forward rus_inverse /usr/local/bin

rus_forward: objects
	$(CC) src/rus_alloc.o src/rus_pars.o src/rus_forward.o $(LFLAGS) -o rus_forward

rus_inverse: objects
	$(CC) src/rus_alloc.o src/rus_pars.o src/xindex.o src/rus_inverse.o $(LFLAGS) -o rus_inverse

uninstall: clean
	-rm -f /usr/local/bin/rus_forward /usr/local/bin/rus_inverse

clean:
	$(MAKE) clean -C src
	-rm -f rus_forward rus_inverse

forward_example: install
	bash example/forward_test.sh

inverse_example: install
	cp example/* /tmp
	rus_inverse

profile: CFLAGS += -pg
profile: clean default

