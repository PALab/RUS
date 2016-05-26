#ifndef RUS_ALLOC_H
#define RUS_ALLOC_H

#include <stddef.h>

void *alloc1(size_t, size_t);
void *realloc1(void*, size_t, size_t);
void **alloc2(size_t, size_t, size_t);
void ***alloc3(size_t, size_t, size_t, size_t);

void free1(void*);
void free2(void**);
void free3(void***);
int *alloc1int(size_t);
int **alloc2int(size_t, size_t);
float *alloc1float(size_t);
float **alloc2float(size_t, size_t);
double *alloc1double(size_t);
double **alloc2double(size_t, size_t);
double ***alloc3double(size_t, size_t, size_t);

void free1int(int*);
void free2int(int**);
void free3int(int***);
void free1float(float*);
void free2float(float**);

void free1double(double*);
void free2double(double**);
void free3double(double***);

/* allocation with error trapping */
void *ealloc1 (size_t n1, size_t size);
void **ealloc2 (size_t n1, size_t n2, size_t size);

int *ealloc1int(size_t n1);
int **ealloc2int(size_t n1, size_t n2);
float *ealloc1float(size_t n1);
float **ealloc2float(size_t n1, size_t n2);

double *ealloc1double(size_t n1);
double **ealloc2double(size_t n1, size_t n2);
double ***ealloc3double(size_t n1, size_t n2, size_t n3);

#endif

