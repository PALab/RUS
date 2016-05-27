#ifndef RUS_PARS_H
#define RUS_PARS_H

#include <stdio.h>

/* GLOBAL DECLARATIONS */
extern int xargc;
extern char **xargv;

/* getpar parameter parsing */
void initargs(int, char**);
int getparint(char*, int*);
int getparuint(char*, unsigned int*);
int getparshort(char*, short*);
int getparushort(char*, unsigned short*);
int getparlong(char*, long*);
int getparulong(char*, unsigned long*);
int getparfloat(char*, float*);
int getpardouble(char*, double*);
int getparstring(char*, char**);
int getparstringarray(char*, char**);
int getnparint(int, char*, int*);
int getnparuint(int, char*, unsigned int*);
int getnparshort(int, char*, short*);
int getnparushort(int, char*, unsigned short*);
int getnparlong(int, char*, long*);
int getnparulong(int, char*, unsigned long *p);
int getnparfloat(int, char*, float*);
int getnpardouble(int, char*, double*);
int getnparstring(int, char*, char**);
int getnparstringarray(int, char*, char**);
int getnpar(int, char*, char*, void*);

/* errors and warnings */
void syserr(char*, ...);
void err(char*, ...);

/* system subroutine calls with error trapping */
FILE *efopen(const char*, const char*);
int efclose(FILE*);
int efseek(FILE*, off_t, int);
off_t eftello(FILE*);
size_t efread(void*, size_t, size_t, FILE*);
size_t efwrite(void*, size_t, size_t, FILE*);

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

void strchop(char*, char*);

/* parameter table */
typedef struct {
	char *name;		/* external name of parameter	*/
	char *asciival;		/* ascii value of parameter	*/
} pointer_table;


#endif

