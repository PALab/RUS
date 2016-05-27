CC = gcc

CFLAGS = -Wall -Werror
LFLAGS = -lm -llapack
FLAGS = $(CFLAGS) $(LFLAGS)

FORWARD = rus_forward
INVERSE = rus_inverse

INCLUDE_PATH = -Iinclude/
SOURCES=$(wildcard src/*.c)
OBJECTS=$(patsubst %.c, %.o, $(SOURCES))
EXECUTABLES = $(FORWARD) $(INVERSE)

default: $(FORWARD) $(INVERSE)

all: install

# forward doesn't require xindex.o
bin/rus_forward: src/rus_alloc.o src/rus_pars.o
	mkdir -p bin
	$(CC) src/rus_forward.c $(INCLUDE_PATH) $(FLAGS) -o bin/rus_forward src/rus_alloc.o src/rus_pars.o

# inverse needs all the objects
bin/rus_inverse: src/rus_alloc.o src/rus_pars.o src/xindex.o
	mkdir -p bin
	$(CC) src/rus_inverse.c $(INCLUDE_PATH) $(FLAGS) -o bin/rus_inverse src/rus_alloc.o src/rus_pars.o src/xindex.o

$(OBJECTS): src/%.o : src/%.c
	$(CC) $(INCLUDE_PATH) $(FLAGS) -c $< -o $@

install: bin/$(FORWARD) bin/$(INVERSE)
	mkdir -p $(HOME)/bin
	cp bin/$(FORWARD) bin/$(INVERSE) $(HOME)/bin

uninstall: clean
	rm -f $(HOME)/bin/$(FORWARD) $(HOME)/bin/$(INVERSE)

clean:
	rm -f $(SRC)/*.o
	rm -f bin/$(FORWARD) bin/$(INVERSE)
	rmdir bin

forward_example: install
	bash example/forward_test.sh

inverse_example: install
	cp example/* /tmp
	$(INVERSE)

