CC = gcc

CFLAGS = -Wall -Werror

LFLAGS = -lm -llapack

TARGET = rus_forward

default: $(TARGET)

all: install

$(TARGET): rus_alloc.o rus_pars.o
	$(CC) $(TARGET).c $(CFLAGS) $(LFLAGS) -o $(TARGET) rus_alloc.o rus_pars.o

rus_alloc.o: rus_alloc.c rus_alloc.h
	$(CC) $(CFLAGS) $(LFLAGS) -c rus_alloc.c

rus_pars.o: rus_pars.c rus_pars.h
	$(CC) $(CFLAGS) $(LFLAGS) -c rus_pars.c

install: $(TARGET)
	mkdir -p $(HOME)/bin
	mv $(TARGET) $(HOME)/bin

uninstall:
	rm -f $(HOME)/bin/$(TARGET)

clean:
	rm -f $(TARGET) *.o
