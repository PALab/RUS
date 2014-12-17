#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h> 

/* Author: Jerome H.L. Le Rousseau (jerome@dix.mines.edu)	        */
/* Center for Wave Phenomena / Physical Acoustic Laboratory             */
/* Colorado School of Mines	 					*/
/* Golden, CO 80401 USA							*/

/* Updated 2014                                         */
/* By: Leighton Watson (lwat054@aucklanduni.ac.nz)      */
/* Physical Acoustic Laboratory                         */
/* University of Auckland, Auckland, 1010, New Zealand  */

/*------------------------ self documentation --------------------------*/
char *sdoc[] = {
"									",
" formod-rp d= d1= d2= d3= rho= lambda= mu= eigenfile=                  ",
"									",
" required parameters:							",
"									",
" d 		order of polynomials used to fit the eigenmodes		",
" d1    	first dimension in cm					",
" d2    	second dimension in cm					",
" d3    	third dimension	in cm					",
" rho 		density in g/cm3					",
"   									",
" ns            number of siffness coefficient given                    ",
"		(determine the symmetry)				",
" hextype       type of hexagonal symmetry: 1=VTI, 2= HTI               ",
"               only required if ns=5                                   ",
" if ns=2	isotropic						",
" c11 									",
" c44									",
" if ns=3	cubic	       						",
" c11 									",
" c12 									",
" c44 									",
" if ns=5 and hextype = 1 VTI       		       			",
" c33									",
" c23 									",
" c12									",
" c44									",
" c66 									",
" if ns=5 and hextype = 2 HTI       		       			",
" c11									",
" c33									",
" c12									",
" c44									",
" c66 									",
" if ns=6   	tetragonal		      				",
" c11									",
" c33									",
" c23 									",
" c12									",
" c44									",
" c66 									",
" if ns=9   	orthorhombic 						",
" c11									",
" c22                                                                   ",
" c33									",
" c23 									",
" c13									",
" c12									",
" c44									",
" c55 									",
" c66 									",
"									",
" optional parameters:							",
"									",
" nfreq=10      gives you the nfreqth first eigen frequnecies 		",
"		on the screen						",
" outeigen=0    if 1 calculates and prints eigenvectors in eigenfile    ", 
" eigenfile= 	file were to put eigenvectors in double format		",
" 									", 
" shape=0      	0 = rectangular parallelepipeds (default)		",
"		1 = ellipsoidal cylinder				",
"		2 = spheroid						",
"									", 
NULL};
/*--------------------- end of self documentation ----------------------*/


/*------------------------ function prototyes --------------------------*/
int delta(int i, int j);

void index_relationship(int *itab, int *ltab, int *mtab, int *ntab,
			int d, int *irk);

double doublefact(int n);

double  volintegral (double d1, double d2, double d3, int l, int m, int n, 
		     int shape);

void e_fill(double *e, int *itab, int *ltab, int *mtab, int *ntab, 
	    int r, double d1, double d2, double d3, double rho, 
	    int shape, int k, int *irk);

void gamma_fill(double *gamma, int *itab, int *ltab, int *mtab, int *ntab, 
		int r, double d1, double  d2, double d3, double ****c, 
		int shape, int k, int *irk);

void stiffness (double ****c, double **cm);

void dqkpart (double a[], int p, int q, int *j, int *k);

void dqkinss (double a[], int p, int q);

void dqksort (int n, double a[]);
/*----------------------------------------------------------------------*/

#define NSTACK 50	/* maximum sort length is 2^NSTACK */
#define NSMALL 7	/* size of array for which insertion sort is fast */
#define FM 7875		/* constants used to generate random pivots */
#define FA 211
#define FC 1663

int main(int argc, char **argv)
{
  /********************* variables declaration **************************/
  int info, itype, lda, ldb, lwork, order; /* variables for lapack function */
  char jobz, uplo; /* variables for lapack function */
  int nfreq; /* number of frequencies displayed on the screen */
  int d; /* dimension of the problem - determine the size r of the partial basis*/
  int shape; /* shape of the body */
  int r; /* actual size of the partial basis */
  int i, j; /* indices */
  int ir1;
  int *itab, *ltab, *mtab, *ntab; /* tabulation of indices */
  int *irk;
  int k;
  int ns; /* symmetry of the system */
  int hextype; /* type of hexagonal symmetry - VTI or HTI*/

  double d1, d2, d3; /* dimension of the sample */
  double rho; /* density */
  double **cm;
  double ****c; /* stiffness tensor */
  double **e, **gamma, *work, **w; /* matrices of the eigenvalue problem */
  double *wsort;
  
  int outeigen; /* 1 if eigenvectors calculated */
  char *eigenfile;

 /** FILE *file; */
  /********************* end variables declaration **********************/
  
  /* hook up getpar to handle the parameters */
  initargs(argc,argv);
  requestdoc(1);
      
  /* get required parameters */
  if (!getparint("d", &d)) err("must specify d!\n");
  if (!getpardouble("d1", &d1)) err("must specify d1!\n");	
  if (!getpardouble("d2", &d2)) err("must specify d2!\n");	
  if (!getpardouble("d3", &d3)) err("must specify d3!\n");	
  if (!getpardouble("rho", &rho)) err("must specify rho!\n");
  if (!getparint("ns", &ns)) err("must specify ns!\n");
  
  cm=ealloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  
  if (ns==2) {
    /* isotropic */
    if (!getpardouble("c11", &cm[0][0])) err("must specify c11!\n");
    if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
    cm[0][0]=cm[0][0]/100;
    cm[3][3]=cm[3][3]/100; 
    cm[1][1]=cm[2][2]=cm[0][0];
    cm[4][4]=cm[5][5]=cm[3][3];	
    cm[0][1]=cm[0][2]=cm[1][2]=cm[0][0]- 2.0*cm[3][3];
    cm[1][0]=cm[2][0]=cm[2][1]=cm[0][0]- 2.0*cm[3][3];

  } else if (ns==3) {
    /* cubic */
    if (!getpardouble("c11", &cm[0][0])) err("must specify c11!\n");
    if (!getpardouble("c12", &cm[0][1])) err("must specify c12!\n");
    if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
    cm[0][0]=cm[0][0]/100;
    cm[0][1]=cm[0][1]/100;
    cm[3][3]=cm[3][3]/100;
    cm[1][1]=cm[2][2]=cm[0][0];	
    cm[4][4]=cm[5][5]=cm[3][3];	
    cm[0][2]=cm[1][2]=cm[0][1];
    cm[2][0]=cm[2][1]=cm[1][0]=cm[0][1];

  } else if (ns==5) {
    /* hexagonal */
    if (!getparint("hextype", &hextype)) err("must specify hextype!\n");

    if (hextype==1) {
      /* VTI */
      if (!getpardouble("c33", &cm[2][2])) err("must specify c33!\n");
      if (!getpardouble("c23", &cm[1][2])) err("must specify c23!\n");
      if (!getpardouble("c12", &cm[0][1])) err("must specify c12!\n");
      if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
      if (!getpardouble("c66", &cm[5][5])) err("must specify c66!\n");

      cm[2][2]=cm[2][2]/100;
      cm[1][2]=cm[1][2]/100;
      cm[0][1]=cm[0][1]/100;
      cm[3][3]=cm[3][3]/100;
      cm[5][5]=cm[5][5]/100;
      cm[0][0]=cm[1][1]=2.0*cm[5][5] + cm[0][1];
      cm[0][2]=cm[2][0]=cm[2][1]=cm[1][2];
      cm[1][0]=cm[0][1];
      cm[4][4]=cm[3][3];

    } else if (hextype==2) {
       
      /* HTI */
      if (!getpardouble("c11", &cm[0][0])) err("must specify c11!\n");
      if (!getpardouble("c33", &cm[2][2])) err("must specify c33!\n");
      if (!getpardouble("c12", &cm[0][1])) err("must specify c12!\n");
      if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
      if (!getpardouble("c66", &cm[5][5])) err("must specify c66!\n");
      cm[0][0]=cm[0][0]/100;
      cm[2][2]=cm[2][2]/100;
      cm[0][1]=cm[0][1]/100;
      cm[3][3]=cm[3][3]/100;
      cm[5][5]=cm[5][5]/100;
      cm[1][2]=cm[2][1]=cm[2][2] - 2.0*cm[3][3];
      cm[0][2]=cm[1][0]=cm[2][0]=cm[0][1];
      cm[1][1]=cm[2][2];
      cm[4][4]=cm[5][5];
      
    }

    else {
      err("for hexagonal symmetry hextype must equal 1 (VTI) or 2 (HTI)!\n");
    }
  }
  
  else if (ns==6){
    /* tetragonal */
    if (!getpardouble("c11", &cm[0][0])) err("must specify c11!\n");
    if (!getpardouble("c33", &cm[2][2])) err("must specify c33!\n");
    if (!getpardouble("c23", &cm[1][2])) err("must specify c23!\n");
    if (!getpardouble("c12", &cm[0][1])) err("must specify c12!\n");
    if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
    if (!getpardouble("c66", &cm[5][5])) err("must specify c66!\n");
    cm[0][0]=cm[0][0]/100;
    cm[2][2]=cm[2][2]/100;
    cm[1][2]=cm[1][2]/100;
    cm[3][3]=cm[3][3]/100;
    cm[0][1]=cm[0][1]/100;
    cm[5][5]=cm[5][5]/100;
    cm[1][1]=cm[0][0];
    cm[0][2]=cm[2][0]=cm[1][2];
    cm[1][0]=cm[0][1];
    cm[2][1]=cm[1][2];
    cm[4][4]=cm[3][3];
  }

  else if (ns==9){/* orthorhombic */
    if (!getpardouble("c11", &cm[0][0])) err("must specify c11!\n");
    if (!getpardouble("c22", &cm[1][1])) err("must specify c22!\n");
    if (!getpardouble("c33", &cm[2][2])) err("must specify c33!\n");
    if (!getpardouble("c23", &cm[1][2])) err("must specify c23!\n");
    if (!getpardouble("c13", &cm[0][2])) err("must specify c13!\n");
    if (!getpardouble("c12", &cm[0][1])) err("must specify c12!\n");
    if (!getpardouble("c44", &cm[3][3])) err("must specify c44!\n");
    if (!getpardouble("c55", &cm[4][4])) err("must specify c55!\n");
    if (!getpardouble("c66", &cm[5][5])) err("must specify c66!\n");
    cm[0][0]=cm[0][0]/100;
    cm[1][1]=cm[1][1]/100;
    cm[2][2]=cm[2][2]/100;
    cm[1][2]=cm[1][2]/100;
    cm[0][2]=cm[0][2]/100;
    cm[0][1]=cm[0][1]/100;
    cm[3][3]=cm[3][3]/100;
    cm[4][4]=cm[4][4]/100;
    cm[5][5]=cm[5][5]/100;
    cm[2][0]=cm[0][2];
    cm[1][0]=cm[0][1];
    cm[2][1]=cm[1][2];
  }

  else err("given elatic moduli does not fit given ns");
  
  

  /* get optional parameters */
  if (!getparint("outeigen", &outeigen)) outeigen=0;
  if (outeigen!=0)
    if (!getparstring("eigenfile", &eigenfile)) 
      err("must specify eigenfile since outeigen>0!\n");
  if (!getparint("shape", &shape)) shape=1; /* changed from zero default to 1 */
  if (!getparint("nfreq", &nfreq)) nfreq=10;

  /* dimension of the problem */
  r= 3*(d+1)*(d+2)*(d+3)/6;
  
  d1=d1/2.0; /* half sample dimensions are used in calculations */
  d2=d2/2.0;
  d3=d3/2.0; 
    
  /* alloc work space*/
  itab=ealloc1int(r);
  ltab=ealloc1int(r);
  mtab=ealloc1int(r);
  ntab=ealloc1int(r);
  
  /* relationship between ir and l,m,n - filling tables */
  irk=ealloc1int(8);
  index_relationship(itab, ltab, mtab, ntab, d, irk); 

  
  
  /* alloc workspace to solve for eigenvalues and eigenfunctions */
  e= (double **) malloc(8*sizeof(double *));
  for (k=0;  k<8; ++k)
    e[k] = ealloc1double(irk[k]*irk[k]);
  
  gamma= (double **) malloc(8*sizeof(double *));
  for (k=0;  k<8; ++k)
    gamma[k] = ealloc1double(irk[k]*irk[k]);
  
  /* filling matrix e */
  for (k=0; k<8; ++k)
    e_fill(e[k], itab, ltab, mtab, ntab, 
	   r, d1, d2, d3, rho, shape, k, irk);
 
  
  /* stiffness tensor calculation*/
  c= (double ****) malloc(sizeof(double ***)*3);
  for (i=0; i<3; ++i)
    c[i]=ealloc3double(3,3,3);
  stiffness (c,  cm);
  
  /* filling matrix gamma  */
  for (k=0; k<8; ++k)
    gamma_fill(gamma[k], itab, ltab, mtab, 
	       ntab, r, d1, d2, d3, c, shape, k, irk);
  

  
  /* clean workspace */
  free1int(itab); 
  free1int(ltab); 
  free1int(mtab); 
  free1int(ntab); 
  for (i=0; i<3; ++i) 
    free3double(c[i]); 
  free(c); 
  fprintf(stderr,"done preparing matrices\n");

  /*-------------------------------------------------------------*/
  /*--------- solve the generalized eigenvalue problem ----------*/
  /*-------------------------------------------------------------*/  
  w= (double **) malloc(sizeof(double *)*8);
  itype=1; 
  if (outeigen==0) jobz='N';
  else jobz='V';
  uplo='U'; 
  for (k=0; k<8; ++k){
    w[k] =ealloc1double(irk[k]);
    lda=ldb=irk[k]; 
    order=irk[k];  
    lwork=MAX(1, 3*order-1);
    work=ealloc1double(lwork);
    /* lapack routine */
    dsygv_(&itype, &jobz, &uplo, &order, gamma[k], 
	   &lda, e[k], &ldb, w[k], work, &lwork, &info);  
    free1double(work);
  } 
  /*-------------------------------------------------------------*/  
  /*-------------------------------------------------------------*/
  /*-------------------------------------------------------------*/
    
  wsort=ealloc1double(r);
   
  for (i=0, k=0; k<8; ++k)
    for (ir1=0;ir1<irk[k];++ir1,++i)
      wsort[i]=w[k][ir1];
   
  /* sorting the eigenfrequencies */
  dqksort(r,wsort);
     
  for (i=0, ir1=0; ir1<nfreq;++i)
    if ((wsort[i]>0) && ((sqrt(wsort[i])/(2.0*PI))>0.00001)){ 
      ++ir1;
      /*fprintf(stderr," f%d = %f\n", ir1, 1000000*sqrt(wsort[i])/(2.0*PI));*/
      fprintf(stderr," f%d = %f\n", ir1, 1000000*sqrt(wsort[i])/(2.0*PI));
      
    }  
  /* modify output of freq values here*/

  
  /* for (k=0;k<8;++k){
    for (ir2=0;ir2<irk[k]*irk[k];++ir2){
      fprintf(stderr,"gamma[%d][%d]=%f\n",k,ir2,gamma[k][ir2]);
        fprintf(stderr,"e[%d][%d]=%f\n",k,ir2,e[k][ir2]);

    }
  }*/     


   /******************* write eigenvectors in files ***************/
  /*if (outeigen==1){
         z=ealloc2double(r,r);  
         for (ir1=0; ir1<r; ++ir1)  
          for (ir2=0; ir2<r; ++ir2)  
	  z[ir2][ir1]=gamma[ir1][ir2*r+ir1];  */
	/* change the order of the array at the same time  */
	/*  since we go from fortran array  */
	/*   to C array */
	/* clean workspace */
	 /*	 free1double(gamma);  
   
        file = efopen(eigenfile, "w"); 
        efwrite(&irf, sizeof(int), 1, file); 
        efwrite(w, sizeof(double), r, file); 
        efwrite(z[0], sizeof(double), r*r, file); 
        efclose(file);*/ 
   /* clean workspace */
    /* free2double(z); */
    /* }*/ 
   
   /* clean workspace */
   /*  free1double(w);  */
   
   /* end of main */
   return EXIT_SUCCESS;
}

int delta(int i, int j){
if (i==j) return 1; 
else return 0;
}



double doublefact(int n){

  if(n==-1) return 1;
  else if (n==0) return 1; 
  else if (n==1) return 1;  
  else return n*doublefact(n-2);
}
double  volintegral (double d1, double d2, double d3, int l, int m, int n, int shape)
{
  if ((l%2==1) || (m%2==1) || (n%2==1)) return 0.0; 
  else 
    switch (shape) {
      /* ell. cylinder shape */
    case 1: return 4.0*PI*pow(d1, l+1)*pow(d2, m+1)*pow(d3, n+1)/(double)(n+1)
	      *doublefact(l-1)*doublefact(m-1)/doublefact(l+m+2);
    /* spheroid shape */
    case 2:  return 4.0*PI*pow(d1, l+1)*pow(d2, m+1)*pow(d3, n+1)
	       *doublefact(l-1)*doublefact(m-1)*doublefact(n-1)/
	       doublefact(l+m+n+3);
    /* rp shape */
    default: return 8.0/((l+1)*(m+1)*(n+1))*pow(d1, l+1)*pow(d2,m+1)*
	       pow(d3, n+1);
    }
}
void e_fill(double *e, int *itab, int *ltab, int *mtab, int *ntab, 
	    int r, double d1, double d2, double d3, double rho, 
	    int shape, int k, int *irk)
{
  int ir1, ir2, i1, i2, l1, l2, m1, m2, n1, n2, irs, irf, ik, irv, irh;
  int l, m, n;
  
  irs=0;
  for (ik=0; ik<k; ++ik)
    irs+=irk[ik];

  irf= irs+irk[k];
    
  for (irv=0, ir1=irs; ir1<irf; ++ir1, ++irv){
    for (irh=0, ir2=irs; ir2<irf; ++ir2, ++irh){
      i1=itab[ir1];
      i2=itab[ir2];
      l1=ltab[ir1];
      l2=ltab[ir2];
      m1=mtab[ir1];
      m2=mtab[ir2];
      n1=ntab[ir1];
      n2=ntab[ir2];

      l=l1+l2;
      m=m1+m2;
      n=n1+n2;
      if (i1!=i2) e[irv*irk[k]+irh]=0.0;
      else e[irv*irk[k]+irh]= rho* volintegral(d1, d2, d3, l, m, n, shape);
    }
  }
}

void stiffness (double ****c, double **cm)
{
  int i, j, k, l;
  int a, b;
 
   for (i=0; i<3; ++i)
    for (j=0; j<3; ++j){
      if ((i==0)&&(j==0)) a=0;
      else if ((i==1)&&(j==1)) a=1;
      else if ((i==2)&&(j==2)) a=2;
      else if (((i==1)&&(j==2))||((i==2)&&(j==1))) a=3;
      else if (((i==0)&&(j==2))||((i==2)&&(j==0))) a=4 ;
      else a=5 ;
      for (k=0; k<3; ++k)
	for (l=0; l<3; ++l){
	  if ((k==0)&&(l==0)) b=0;
	  else if ((k==1)&&(l==1)) b=1;
	  else if ((k==2)&&(l==2)) b=2;
	  else if (((k==1)&&(l==2))||((k==2)&&(l==1))) b=3;
	  else if (((k==0)&&(l==2))||((k==2)&&(l==0))) b=4 ;
	  else b=5 ;
	  c[i][j][k][l]=cm[a][b];
 	}
    }
 } 

void gamma_fill(double *gamma, int *itab, int *ltab, int *mtab, int *ntab, 
		int r, double d1, double  d2, double d3, double ****c, 
		int shape, int k, int *irk)
{
  int ir1, ir2, i1, i2, l1, l2, m1, m2, n1, n2, irs, irf, ik, irv, irh;
  int j1, j2;
  int l, m, n;
  irs=0;
  
  for (ik=0; ik<k; ++ik)
    irs+=irk[ik];

  irf= irs+irk[k];

  for (irv=0, ir1=irs; ir1<irf; ++ir1, ++irv){
    for (irh=0, ir2=irs; ir2<irf; ++ir2, ++irh){
      i1=itab[ir1];
      i2=itab[ir2];
      l1=ltab[ir1];
      l2=ltab[ir2];
      m1=mtab[ir1];
      m2=mtab[ir2];
      n1=ntab[ir1];
      n2=ntab[ir2];
      /* initialize the value of gamma to 0 */
      gamma[irv*irk[k]+irh]=0.0;
      if (l1>0) {
	j1=0;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-2;
	  m=m1+m2;
	  n=n1+n2;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*l1*l2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2-1;
	  m=m1+m2-1;
	  n=n1+n2;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*l1*m2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2-1;
	  m=m1+m2;
	  n=n1+n2-1;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*l1*n2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
      }
      if (m1>0) {
	j1=1;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-1;
	  m=m1+m2-1;
	  n=n1+n2;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*m1*l2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2;
	  m=m1+m2-2;
	  n=n1+n2; 
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*m1*m2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2;
	  m=m1+m2-1;
	  n=n1+n2-1;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*m1*n2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
      }
      if (n1>0) {
	j1=2;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-1;
	  m=m1+m2;
	  n=n1+n2-1;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*n1*l2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2;
	  m=m1+m2-1;
	  n=n1+n2-1;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*n1*m2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2;
	  m=m1+m2;
	  n=n1+n2-2;
	  gamma[irv*irk[k]+irh] +=c[i1][j1][i2][j2]*n1*n2*
	    volintegral(d1, d2, d3, l, m, n, shape);
	}
      }
    }
  }
}



void dqkpart (double a[], int p, int q, int *j, int *k)
{
	int pivot,left,right;
	double apivot,temp;
	static long int seed=0L;
 
	/* choose random pivot element between p and q, inclusive */
	seed = (seed*FA+FC)%FM;
	pivot = p+(q-p)*(double)seed/(double)FM;
	if (pivot<p) pivot = p;
	if (pivot>q) pivot = q;
	apivot = a[pivot];

	/* initialize left and right pointers and loop until break */
	for (left=p,right=q;;) {
		/*
		 * increment left pointer until either
		 * (1) an element greater than the pivot element is found, or
		 * (2) the upper bound of the input subarray is reached
		 */
		while (a[left]<=apivot && left<q) left++;
 
		/*
		 * decrement right pointer until either
		 * (1) an element less than the pivot element is found, or
		 * (2) the lower bound of the input subarray is reached
		 */
		while (a[right]>=apivot && right>p) right--;
 
		/* if left pointer is still to the left of right pointer */
		if (left<right) {
			/* exchange left and right elements */
			temp = a[left];
			a[left++] = a[right];
			a[right--] = temp;
		} 
		/* else, if pointers are equal or have crossed, break */
		else break;
	}
	/* if left pointer has not crossed pivot */
	if (left<pivot) {

		/* exchange elements at left and pivot */
		temp = a[left];
		a[left++] = a[pivot];
		a[pivot] = temp;
	}
	/* else, if right pointer has not crossed pivot */
	else if (pivot<right) {

		/* exchange elements at pivot and right */
		temp = a[right];
		a[right--] = a[pivot];
		a[pivot] = temp;
	}
	/* left and right pointers have now crossed; set output bounds */
	*j = right;
	*k = left;
}

void dqkinss (double a[], int p, int q)
{
	int i,j;
	double ai;

	for (i=p+1; i<=q; i++) {
		for (ai=a[i],j=i; j>p && a[j-1]>ai; j--)
			a[j] = a[j-1];
		a[j] = ai;
	}
}

void dqksort (int n, double a[]){
	int pstack[NSTACK],qstack[NSTACK],j,k,p,q,top=0;

	/* initialize subarray lower and upper bounds to entire array */
	pstack[top] = 0;
	qstack[top++] = n-1;

	/* while subarrays remain to be sorted */
	while(top!=0) {

		/* get a subarray off the stack */
		p = pstack[--top];
		q = qstack[top];

		/* while subarray can be partitioned efficiently */
		while(q-p>NSMALL) {

			/* partition subarray into two subarrays */
			dqkpart(a,p,q,&j,&k);

			/* save larger of the two subarrays on stack */
			if (j-p<q-k) {
				pstack[top] = k;
				qstack[top++] = q;
				q = j;
			} else {
				pstack[top] = p;
				qstack[top++] = j;
				p = k;
			}
		}
		/* use insertion sort to finish sorting small subarray */
		dqkinss(a,p,q);
	}
}
void index_relationship(int *itab, int *ltab, int *mtab, int *ntab,
			int d, int *irk)
{
  int ir, i, l, m, n;
  ir=0;
  
  /* k=0 */
  irk[0]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[0];
		}
	  if (i==1)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[0];
		} 
	  if (i==2)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[0];
		}
	}
  fprintf(stderr, "irk[0]=%d\n", irk[0]);
  
  /* k=1 */
  irk[1]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[1];
		}
	  if (i==1)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[1];
		} 
	  if (i==2)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[1];
		}
	}
  fprintf(stderr, "irk[1]=%d\n", irk[1]);
  
  /* k=2 */
  irk[2]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[2];
		}
	  if (i==1)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[2];
		} 
	  if (i==2)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[2];
		}
	}
  fprintf(stderr, "irk[2]=%d\n", irk[2]);
  
  /* k=3 */
  irk[3]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[3];
		}
	  if (i==1)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[3];
		} 
	  if (i==2)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[3];
		}
	}
  fprintf(stderr, "irk[3]=%d\n", irk[3]);
  
  /* k=4 */
  irk[4]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[4];
		}
	  if (i==1)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[4];
		}
		 
	  if (i==2)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[4];
		}	
	}
  fprintf(stderr, "irk[4]=%d\n", irk[4]);
  
 /* k=5 */
  irk[5]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==1)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[5];
		}
	  
	  if (i==1)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[5];
		} 
	  if (i==2)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[5];
		}
	}
  fprintf(stderr, "irk[5]=%d\n", irk[5]);
  
   /* k=6 */
  irk[6]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[6];
		}
	  if (i==1)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[6];
		} 
	  if (i==2)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[6];
		}
	}
  fprintf(stderr, "irk[6]=%d\n", irk[6]);
  
   /* k=7 */
  irk[7]=0;
  for (i=0; i<3; ++i)
    for (l=0; l<=d; ++l)
      for (m=0; m+l<=d; ++m)
	for (n=0; n+l+m<=d; ++n){
	  if (i==0)
	    if (l%2==1)
	      if (m%2==1)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;

		  ntab[ir]=n;
		  ++ir;
		  ++irk[7];
		}
	  if (i==1)
	    if (l%2==0)
	      if (m%2==0)
		if(n%2==1){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[7];
		} 
	  if (i==2)
	    if (l%2==0)
	      if (m%2==1)
		if(n%2==0){
		  itab[ir]=i; 
		  ltab[ir]=l;
		  mtab[ir]=m;
		  ntab[ir]=n;
		  ++ir;
		  ++irk[7];
		}
	}
  fprintf(stderr, "irk[7]=%d\n", irk[7]);
  
}
