#include "rus_pars.h"
#include "rus_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>

#define ERROR NULL

#ifndef FALSE
#define FALSE (0)
#endif

#ifndef TRUE
#define TRUE (1)
#endif

#define PAR_NAMES_MAX 512

/* GLOBAL DECLARATIONS */
extern int xargc;
extern char **xargv;

/* global variables declared and used internally */
static pointer_table *argtbl;	/* parameter table		*/
static int nargs;		/* number of args that parse	*/
static int tabled = FALSE;	/* true when parameters tabled 	*/
static size_t targc;		/* total number of args		*/
static char **targv;		/* pointer to arg strings	*/
static char *argstr;		/* storage for command line	*/

/* functions declared and used internally */
static int getparindex (int n, char *name);
static void getparinit(void);
static void tabulate (size_t argc, char **argv);
static char *getpfname (void);

static char* par_names[PAR_NAMES_MAX];
static int par_count=0;
static int parcheck = 0;

/* make command line args available to subroutines -- re-entrant version */
void initargs(int argc, char **argv)
{
	memset( par_names ,0 ,sizeof(par_names) );
	par_names[0] = "par";
	par_names[1] = "lheader";
	par_count=2;

	xargc = argc; xargv = argv;
	if(tabled==TRUE){
		free(argstr);
		free(targv);
		free(argtbl);
	}
	tabled =  FALSE;
	return;
}

size_t efread(void *bufptr, size_t size, size_t count, FILE *stream)
{
	size_t nread;

	if (!size) err("%s: efread: fread given 0 item size", __FILE__);

	nread = fread(bufptr, size, count, stream);

	if (nread != count && ferror(stream))
		      err("%s: efread: fread only %d items of %d",
				__FILE__, nread, count);

	return nread;
}

FILE *efopen(const char *file, const char *mode)
{
	FILE *stream;

	if (ERROR == (stream = fopen(file, mode)))
		err("%s: efopen: fopen failed", __FILE__);
	
	return stream;
}


int efclose(FILE *stream)
{
	int status;

	if (EOF == (status = fclose(stream)))
		      err("%s: efclose: fclose failed", __FILE__);

	return status;
}


int efseek(FILE *stream, off_t offset, int origin)
{
	if (fseek(stream, offset, origin))  /* non-zero => error */
		      err("%s: efseek: fseek failed", __FILE__);

	return 0;
}



off_t eftello(FILE *streem)
{
	off_t eposition;
	off_t test=-1;

	eposition = ftello(streem);
	if (test == eposition) {
		fprintf(stderr,"sizeof(off_t)=%lu\n",
				(unsigned long) sizeof(eposition));
	}
	

	return eposition;
}

/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} */
short eatoh(char *s)
{
	long n = strtol(s, NULL, 10);
	
	if ( (n > SHRT_MAX) || (n < SHRT_MIN) || (errno == ERANGE) )
		err("%s: eatoh: overflow", __FILE__);

	return (short) n;
}


/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} */
unsigned short eatou(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if ( (n > USHRT_MAX) || (errno == ERANGE) )
		err("%s: eatou: overflow", __FILE__);

	return (unsigned short) n;
}


/* eatoi - convert string s to integer {INT_MIN:INT_MAX} */
int eatoi(char *s)
{
	long n = strtol(s, NULL, 10);

	if ( (n > INT_MAX) || (n < INT_MIN) || (errno == ERANGE) )
		err("%s: eatoi: overflow", __FILE__);

	return (int) n;
}


/* eatop - convert string s to unsigned integer {0:UINT_MAX} */
unsigned int eatop(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if ( (n > UINT_MAX) || (errno == ERANGE) )
		err("%s: eatop: overflow", __FILE__);

	return (unsigned int) n;
}


/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX} */
long eatol(char *s)
{
	long n = strtol(s, NULL, 10);

	if (errno == ERANGE)
		err("%s: eatol: overflow", __FILE__);

	return n;
}


/* eatov - convert string s to unsigned long {0:ULONG_MAX} */
unsigned long eatov(char *s)
{
	unsigned long n = strtoul(s, NULL, 10);

	if (errno == ERANGE)
		err("%s: eatov: overflow", __FILE__);

	return n;
}


/* eatof - convert string s to float {-FLT_MAX:FLT_MAX} */
float eatof(char *s)
{
	float x = strtod(s, NULL);

	if ( (x > FLT_MAX) || (x < -FLT_MAX) || (errno == ERANGE) )
		err("%s: eatof: overflow", __FILE__);

	return (float) x;
}


/* eatod - convert string s to double {-DBL_MAX:DBL_MAX} */
double eatod(char *s)
{
	double x = strtod(s, NULL);

	/* errno == ERANGE suffices if compiler sets errno on overflow */
	if ( (errno == ERANGE) || (x > DBL_MAX) || (x < -DBL_MAX) )
		err("%s: eatod: overflow", __FILE__);

	return x;
}

/* functions to get values for the last occurrence of a parameter name */
int getparint (char *name, int *ptr)
{
	return getnpar(0,name,"i",ptr);
}
int getparuint (char *name, unsigned int *ptr)
{
	return getnpar(0,name,"p",ptr);
}
int getparshort (char *name, short *ptr)
{
	return getnpar(0,name,"h",ptr);
}
int getparushort (char *name, unsigned short *ptr)
{
	return getnpar(0,name,"u",ptr);
}
int getparlong (char *name, long *ptr)
{
	return getnpar(0,name,"l",ptr);
}
int getparulong (char *name, unsigned long *ptr)
{
	return getnpar(0,name,"v",ptr);
}
int getparfloat (char *name, float *ptr)
{
	return getnpar(0,name,"f",ptr);
}
int getpardouble (char *name, double *ptr)
{
	return getnpar(0,name,"d",ptr);
}
int getparstring (char *name, char **ptr)
{
	return getnpar(0,name,"s",ptr);
}
int getparstringarray (char *name, char **ptr)
{
	return getnpar(0,name,"a",ptr);
}
int getpar (char *name, char *type, void *ptr)
{
	return getnpar(0,name,type,ptr);
}

/* functions to get values for the n'th occurrence of a parameter name */
int getnparint (int n, char *name, int *ptr)
{
	return getnpar(n,name,"i",ptr);
}
int getnparuint (int n, char *name, unsigned int *ptr)
{
	return getnpar(n,name,"p",ptr);
}
int getnparshort (int n, char *name, short *ptr)
{
	return getnpar(n,name,"h",ptr);
}
int getnparushort (int n, char *name, unsigned short *ptr)
{
	return getnpar(n,name,"u",ptr);
}
int getnparlong (int n, char *name, long *ptr)
{
	return getnpar(n,name,"l",ptr);
}
int getnparulong (int n, char *name, unsigned long *ptr)
{
	return getnpar(n,name,"v",ptr);
}
int getnparfloat (int n, char *name, float *ptr)
{
	return getnpar(n,name,"f",ptr);
}
int getnpardouble (int n, char *name, double *ptr)
{
	return getnpar(n,name,"d",ptr);
}
int getnparstring (int n, char *name, char **ptr)
{
	return getnpar(n,name,"s",ptr);
}
int getnparstringarray (int n, char *name, char **ptr)
{
	return getnpar(n,name,"a",ptr);
}
int getnpar (int n, char *name, char *type, void *ptr)
{
	int i;			/* index of name in symbol table	*/
	int j;		  /* index for par_names[]		*/
	int nval;		/* number of parameter values found	*/
	char *aval;		/* ascii field of symbol		*/

/*--------------------------------------------------------------------*\
   getpar gets called in loops reading traces in some programs.  So
   check for having seen this name before. Also make sure we don't
   walk off the end of the table.
\*--------------------------------------------------------------------*/

	if( parcheck && strcmp( "lheader" ,name ) ){
	   fprintf( stderr ,"getpar() call after checkpars(): %s\n" ,name );
	}

	for( j=0; j<par_count; j++ ){
	   if( !strcmp( par_names[j] ,name ) ){
		break;
	   }
	}

	if( j >= par_count && par_count < PAR_NAMES_MAX ){
	   par_names[par_count++] = name;
	}

	if(  par_count == PAR_NAMES_MAX ){
	   fprintf( stderr, " %s exceeded PAR_NAMES_MAX %d \n" ,xargv[0] ,PAR_NAMES_MAX );
	}

	if (xargc == 1) return 0;
	if (!tabled) getparinit();/* Tabulate command line and parfile */
	i = getparindex(n,name);/* Get parameter index */
	if (i < 0) return 0;	/* Not there */

	if (0 == ptr) {
	   err("%s: getnpar called with 0 pointer, type = %s", __FILE__,type);
	}
	  

	/*
	 * handle string type as a special case, since a string
	 * may contain commas.
	 */
	if (type[0]=='s') {
		*((char**)ptr) = argtbl[i].asciival;
		return 1;
	}

	/* convert vector of ascii values to numeric values */
	for (nval=0,aval=argtbl[i].asciival; *aval; nval++) {
		switch (type[0]) {
			case 'i':
				*(int*)ptr = eatoi(aval);
				ptr = (int*)ptr+1;
				break;
			case 'p':
				*(unsigned int*)ptr = eatop(aval);
				ptr = (unsigned int*)ptr+1;
				break;
			case 'h':
				*(short*)ptr = eatoh(aval);
				ptr = (short*)ptr+1;
				break;
			case 'u':
				*(unsigned short*)ptr = eatou(aval);
				ptr = (unsigned short*)ptr+1;
				break;
			case 'l':
				*(long*)ptr = eatol(aval);
				ptr = (long*)ptr+1;
				break;
			case 'v':
				*(unsigned long*)ptr = eatov(aval);
				ptr = (unsigned long*)ptr+1;
				break;
			case 'f':
				*(float*)ptr = eatof(aval);
				ptr = (float*)ptr+1;
				break;
			case 'd':
				*(double*)ptr = eatod(aval);
				ptr = (double*)ptr+1;
				break;
			case 'a':
				{ char *tmpstr="";
				   tmpstr = ealloc1(strlen(aval)+1,1);

				   strchop(aval,tmpstr);
				   *(char**)ptr = tmpstr;
				   ptr=(char **)ptr + 1;
				}
				   break;
			default:
				err("%s: invalid parameter type = %s",
					__FILE__,type);
		}
		while (*aval++ != ',') {
			if (!*aval) break;
		}
	}
	return nval;
}

/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex (int n, char *name)
{
	int i;
	if (n==0) {
		for (i=nargs-1; i>=0; --i)
			if (!strcmp(name,argtbl[i].name)) break;
		return i;
	} else {
		for (i=0; i<nargs; ++i)
			if (!strcmp(name,argtbl[i].name))
				if (--n==0) break;
		if (i<nargs)
			return i;
		else
			return -1;
	}
}

/* Initialize getpar */

static void getparinit (void)
{
	static char *pfname;	/* name of parameter file		*/
	FILE *pffd=NULL;	/* file id of parameter file		*/
	size_t pflen;		/* length of parameter file in bytes	*/
	/*static size_t pfargc;*/	/* arg count from parameter file	*/
	int parfile;		/* parfile existence flag		*/
	int argstrlen=0;
	char *pargstr;		/* storage for parameter file args	*/
	size_t nread=0;		/* bytes fread				*/
	int i, j;		/* counters				*/
	int start = TRUE;
	int debug = FALSE;
	int quote = FALSE;


	tabled = TRUE;		/* remember table is built		*/


	/* Check if xargc was initiated */

	if(!xargc)
		err("%s: xargc=%d -- not initiated in main", __FILE__, xargc);

	/* Space needed for command lines */

	for (i = 1, argstrlen = 0; i < xargc; i++) {
		argstrlen += strlen(xargv[i]) + 1;
	}

	/* Get parfile name if there is one */

	if ((pfname = getpfname())) {
		parfile = TRUE;
	} else {
		parfile = FALSE;
	}

	if (parfile) {
	 	pffd = efopen(pfname, "r");

		/* Get the length */
		efseek(pffd, 0, SEEK_END);

		pflen = (off_t) eftello(pffd);

		rewind(pffd);
		argstrlen += pflen;
	} else {
		pflen = 0;
	}

/*--------------------------------------------------------------------*\
   Allocate space for command line and parameter file. The pointer
   table could be as large as the string buffer, but no larger.

   The parser logic has been completely rewritten to prevent bad
   input from crashing the program.

   Reginald H. Beardsley			    rhb@acm.org
\*--------------------------------------------------------------------*/

	argstr = (char *) ealloc1(argstrlen+1, 1);
	targv = (char **) ealloc1((argstrlen+1)/4,sizeof(char*));

	if (parfile) {
		/* Read the parfile */

		nread = efread(argstr, 1, pflen, pffd);
  		if (nread != pflen) {
  	 	    err("%s: fread only %d bytes out of %d from %s",
  					__FILE__,  nread, pflen, pfname);
		}
		efclose(pffd);


	} else {
		/* pfargc = 0; */
	}


	/* force input to valid 7 bit ASCII */

	for( i=0; i<nread; i++ ){
	    argstr[i] &= 0x7F;
	}

	/* tokenize the input */

	j = 0;

	for( i=0; i<nread; i++ ){

	    /* look for start of token */

	    if( start ){

 /* getpars.c:475: warning: subscript has type `char' */
		if( isgraph( (int)argstr[i] ) ){
		    targv[j] = &(argstr[i]);
		    start = !start;
		    j++;

		}else{
		    argstr[i] = 0;

		}

	    /* terminate token */

/* getpars.c:487: warning: subscript has type `char' */
	    }else if( !quote && isspace( (int)argstr[i] ) ){
		argstr[i] = 0;
		start = !start;

	    }

	    /* toggle quote semaphore */

	    if( argstr[i] == '\'' || argstr[i] == '\"' ){
		quote = !quote;

	    }

	}

	/* display all tokens */

	if( debug ){

	    i=0;
	    while( i < j && targv[i] != 0 ){
		if( strlen( targv[i] ) ){
		    fprintf( stderr ,"%d -> %s\n" ,i ,targv[i] );
		}
		i++;

	    }
	}

	/* discard non-parameter tokens */

	i=0;
	targc=0;
	while( i < j && targv[i] != 0 ){
	    if( strchr( targv[i] ,'=' ) ){
		targv[targc] = targv[i];
		targc++;
	    }
	    i++;
	}

	/* Copy command line arguments */

	for (j = 1, pargstr = argstr + pflen + 1; j < xargc; j++) {
		strcpy(pargstr,xargv[j]);
		targv[targc++] = pargstr;
		pargstr += strlen(xargv[j]) + 1;
	}

	/* Allocate space for the pointer table */

	argtbl = (pointer_table*) ealloc1(targc, sizeof(pointer_table));

	/* Tabulate targv */

	tabulate(targc, targv);

	return;
}
#define PFNAME "par="
/* Get name of parameter file */
static char *getpfname (void)
{
	int i;
	size_t pfnamelen;

	pfnamelen = strlen(PFNAME);
	for (i = xargc-1 ; i > 0 ; i--) {
		if(!strncmp(PFNAME, xargv[i], pfnamelen)
		    && strlen(xargv[i]) != pfnamelen) {
			return xargv[i] + pfnamelen;
		}
	}
	return NULL;
}

#define iswhite(c)	((c) == ' ' || (c) == '\t' || (c) == '\n')

/* Install symbol table */
static void tabulate (size_t argc, char **argv)
{
	int i;
	char *eqptr;
	int debug=FALSE;

	for (i = 0, nargs = 0 ; i < argc; i++) {
		eqptr = strchr(argv[i], '=');
		if (eqptr) {
			argtbl[nargs].name = argv[i];
			argtbl[nargs].asciival = eqptr + 1;
			*eqptr = (char)0;

			/* Debugging dump */
			if( debug ){
				fprintf(stderr,
				"argtbl[%d]: name=%s asciival=%s\n",
				nargs,argtbl[nargs].name,argtbl[nargs].asciival);

			}
			nargs++;
		}
	}
	return;
}

void syserr(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nsyserr: fflush failed on stdout");
	}
	fprintf(stderr, "\n%s: ", xargv[0]);
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, " (%s)\n", strerror(errno));
	exit(EXIT_FAILURE);
}

void strchop(char *s, char *t)
/***********************************************************************
strchop - chop off the tail end of a string "s" after a "," returning
	  the front part of "s" as "t".
************************************************************************
Notes:
Based on strcpy in Kernighan and Ritchie's C [ANSI C] book, p. 106.
************************************************************************
Author: CWP: John Stockwell and Jack K. Cohen, July 1995
***********************************************************************/
{

	while ( (*s != ',') && (*s != '\0') ) {
		 *t++ = *s++;
	}
	*t='\0';
}

