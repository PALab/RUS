#include "rus_alloc.h"
#include "rus_pars.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#define ERROR NULL

/* allocate a 1-d array */
void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	size_t i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}

	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
int **alloc2int(size_t n1, size_t n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(size_t n1)
{
	return (double*)alloc1(n1,sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, size_t n1)
{
	return (double*)realloc1(v,n1,sizeof(double));
}


/* free a 1-d array of doubles */
void free1double(double *p)
{
	free1(p);
}

/* allocate a 2-d array of doubles */
double **alloc2double(size_t n1, size_t n2)
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(size_t n1, size_t n2, size_t n3)
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array */
void *ealloc1 (size_t n1, size_t size)
{
	void *p;

	if (ERROR == (p=alloc1(n1, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 2-d array */
void **ealloc2 (size_t n1, size_t n2, size_t size)
{
	void **p;

	if (ERROR == (p=alloc2(n1, n2, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 1-d array of ints */
int *ealloc1int(size_t n1)
{
	int *p;

	if (ERROR == (p=alloc1int(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 2-d array of ints */
int **ealloc2int(size_t n1, size_t n2)
{
	int **p;

	if (ERROR == (p=alloc2int(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 1-d array of floats */
float *ealloc1float(size_t n1)
{
	float *p;

	if (ERROR == (p=alloc1float(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 2-d array of floats */
float **ealloc2float(size_t n1, size_t n2)
{
	float **p;

	if (ERROR == (p=alloc2float(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 1-d array of doubles */
double *ealloc1double(size_t n1)
{
	double *p;

	if (ERROR == (p=alloc1double(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 2-d array of doubles */
double **ealloc2double(size_t n1, size_t n2)
{
	double **p;

	if (ERROR == (p=alloc2double(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 3-d array of doubles */
double ***ealloc3double(size_t n1, size_t n2, size_t n3)
{
	double ***p;

	if (ERROR == (p=alloc3double(n1, n2, n3)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

