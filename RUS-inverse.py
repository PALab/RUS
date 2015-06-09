# Author: Jerome H.L. Le Rousseau (jerome@dix.mines.edu)
# Center for Wave Phenomena / Physical Acoustic Laboratory
# Colorado School of Mines
# Golden, CO 80401 USA

# Updated 2014
# By: Leighton Watson (lwat054@aucklanduni.ac.nz)
# Physical Acoustic Laboratory
# University of Auckland, Auckland, 1010, New Zealand

# Translated to Python 2015
# By: Paul Freeman (pfre484@aucklanduni.ac.nz)
# Computer Science Department
# University of Auckland, Auckland, 1010, New Zealand

import argparse
import sys
import scipy

NSTACK = 50	# maximum sort length is 2^NSTACK
NSMALL = 7	# size of array for which insertion sort is fast
FM = 7875	# constants used to generate random pivots
FA = 211
FC = 1663

parser = argparse.ArgumentParser(description='Inverse Algorithm')
parser.add_argument(
	'--d',
	type=int,
	required=True,
	help='order of polynomials used to estimate the eigenvectors')
parser.add_argument(
	'--shape',
	nargs='?',
	const='1',
	default='1',
	type=int,
	choices=[0,1,2],
	help='0=sphere, 1=cylinder, 2=parallelepiped')
parser.add_argument(
	'--ns',
	type=int,
	required=True,
	choices=[2,3,5,6,9],
	help='number of cijs')
parser.add_argument(
	'--hextype',
	nargs=1,
	type=int,
	choices=[1,2],
	help='hextype - 1=VTI, 2=HTI. Type of hexagonal symetry (Only matters for ns=5)')
parser.add_argument(
	'--d1f',
	type=float,
	required=True,
	help='dimension 1 in cm (diameter for cyl. or sphere)')
parser.add_argument(
	'--d2f',
	type=float,
	required=True,
	help='dimension 2 in cm (diameter for cyl. or sphere)')
parser.add_argument(
	'--d3f',
	type=float,
	required=True,
	help='dimension 3 in cm (height for cyl. diameter for sphere)')
parser.add_argument(
	'--rhof',
	type=float,
	required=True,
	help='density in grams/cm^3')
parser.add_argument(
	'--freqmin',
	type=float,
	required=True
	help='lower frequency bound for inversion in MHz (set >1 KHz lower than your lowest measured value)')
parser.add_argument(
	'--freqmax',
	type=float,
	required=True
	help='upper frequency bound for inversion in MHz (set >5 or 10KHz higher than your highest value used for THIS particular fit as defined by Line 1 of freq_data)')
#parser.add_argument('--c11', nargs=1, type=float)
#parser.add_argument('--c12', nargs=1, type=float)
#parser.add_argument('--c13', nargs=1, type=float)
#parser.add_argument('--c22', nargs=1, type=float)
#parser.add_argument('--c23', nargs=1, type=float)
#parser.add_argument('--c33', nargs=1, type=float)
#parser.add_argument('--c44', nargs=1, type=float)
#parser.add_argument('--c55', nargs=1, type=float)
#parser.add_argument('--c66', nargs=1, type=float)

args = parser.parse_args()

# matrices of the eigenvalue problem and for the function gradiant
measurement = "freq_data"  # CHANGE THIS LINE TO APPROPRIATE DIRECTORY
     
guessf = scipy.zeros((args.ns))

  /* print out initial guess */
  for (is=0;is<ns;++is){
    fscanf(parameterfile, "%f",&guessf[is]); 
   /*  fprintf(stderr, "GUESS #1=%f GPa\n",guessf[is]);*/
    fprintf(stderr, "%f\n",guessf[is]);
    guessf[is]=guessf[is]/100;
  }
  fclose (parameterfile);
  
  d1=d1f;
  d2=d2f;
  d3=d3f;
  rho=rhof;
  
  /* dimension of the problem */
  r= 3*(d+1)*(d+2)*(d+3)/6;

  d1=d1/2.0; /* half sample dimensions are used in calculations */
  d2=d2/2.0;
  d3=d3/2.0; 
  
  guess=ealloc1double(ns);
  for (is=0;is<ns;++is)
    guess[is]=guessf[is];
  free1float(guessf);

  /* get measured frequencies from file */
  measurementfile = fopen(measurement, "r");
  fscanf(measurementfile, "%d", &nfreq);
  fprintf(stderr, "nfreq=%d\n ", nfreq);
  freq=alloc1float(nfreq);
  weight=alloc1float(nfreq);
  for(ifreq=0; ifreq<nfreq; ++ifreq){
    fscanf(measurementfile, "%f", &freq[ifreq]);
    /* scale freqs back to MHz from Hz*/
    /* freq[ifreq]=freq[ifreq]/1000; */ 
    fprintf(stderr, "freq=%f\n ", freq[ifreq]);
    fscanf(measurementfile, "%f", &weight[ifreq]);
  }
  fclose(measurementfile); 

  
  


  /* alloc work space*/
  itab=alloc1int(r);
  ltab=alloc1int(r);
  mtab=alloc1int(r);
  ntab=alloc1int(r);

  /* relationship between ir and l,m,n - filling tables */
  irk=alloc1int(8);
  index_relationship(itab, ltab, mtab, ntab, d, irk); 
  /* ifr1 = 0 and ifr2 = ndata, number of frequencies used for inversion */
  xindex(nfreq, freq, freqmin, &ifr1);
  xindex(nfreq, freq, freqmax, &ifr2);
  ifr2++; /* add one to ifr2 so the subtraction ifr2 - ifr1 gives the number of frequencies ndata*/
  ndata=ifr2-ifr1; /* ndata = number of frequencies used for inversion - first line from freq_data */
  fprintf(stderr,"ndata=%d\n",ndata);
  y=ealloc1double(ndata);
  for (i=0, ifr=ifr1; ifr<ifr2; ++ifr, ++i){
    y[i]=freq[ifr]; /* set y equal to freq_data values. Therefore, y is an array of length ndata */
  }
  sig=ealloc1double(ndata);
  for (i=0, ifr=ifr1; ifr<ifr2; ++ifr, ++i)
    sig[i]=weight[ifr]; /* set sig equal to the weightings from the second column of freq_data. Is an array of length ndata */
  /* sig is not a true 'error' but the weighting from freq_data */
  /* is multipled by the difference in the misfit function/chisq so modes with a zero weighting are ignored */
  
  
  ia=ealloc1int(ndata);
  for (i=0;i<ndata;++i)
    ia[i]=1.0; /* creates a array of length ndata with entries of 1 */
  /* was having issues with this before. The issues were with determining the size of ia due to incorrect use sizeof*/
  /* ia is specified as a pointer. Cannot determine length of a pointer to an array */
  
  covar=ealloc2double(ns,ns);
  for (is1=0;is1<ns;++is1)
    for (is2=0;is2<ns;++is2)
      covar[is1][is2]=delta(is1,is2); /* identity matrix. Size of ns/number of cijs */
     
  alpha=ealloc2double(ns,ns);
  for (is1=0;is1<ns;++is1)
    for (is2=0;is2<ns;++is2)
      alpha[is1][is2]=delta(is1,is2); /* identity matrix. Size of ns/number of cijs */
  chisq=0.0;
  alamda=-1.0;

  niter=100; /* set number of iterations */
  for (iter=0;iter<niter;++iter){
    mrqmin(d,r,itab,ltab,mtab,ntab,
	   irk,d1,d2,d3,rho,shape,freqmin,
	   y,sig,ndata,guess,ia,ns,covar,alpha,&chisq, hextype, formod,&alamda);
    fprintf(stderr,"iter #%i\n",iter); /* print iteration number */
    for (is=0; is<ns;++is) /* ns = dimension of symmetry */
      /*fprintf(stderr,"guess=%f\n",100*guess[is]); */
      fprintf(stderr,"%f\n",100*guess[is]); /* print estimated cij values*/
    
  }
    
 /* end of main */

  return EXIT_SUCCESS;
}

void formod(int d,int r,int *itab,int *ltab,int *mtab,int *ntab, 
	    int *irk,double d1,double d2,double d3, 
	    double rho,int shape,float freqmin,int ndata,
	    double *a, double *y, double **dyda, int ns, int hextype)
{
int info, itype, lda, ldb, lwork, order; /* variables for lapack function */
  char jobz, uplo; /* variables for lapack function */
  int i,j; /* indices */
  int ir1, ir2, irf;
  int ifw, iw; 
  int k;

  double **cm;
  double ****c; /* stiffness tensor and derivatives*/
  double **e, **gamma, *work, **w,  **z; 
  double *wsort;
  double *wnosort;
  int *indice;
  char freqs[]="/home/paul/RUS/example/predictedf"; /*CHANGE THIS LINE TO APPROPRIATE DIRECTORY*/
  FILE *freqfile;
  
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  if (ns==2) {
    /* isotropic */
    
    cm[0][0]=a[0];
    cm[3][3]=a[1];
    
    cm[1][1]=cm[2][2]=cm[0][0];
    cm[4][4]=cm[5][5]=cm[3][3];	
    cm[0][1]=cm[0][2]=cm[1][2]=cm[0][0]-2.0*cm[3][3];
    cm[1][0]=cm[2][0]=cm[2][1]=cm[0][0]-2.0*cm[3][3];

  } else if (ns==3) {
    /* cubic */
	
    cm[0][0]=a[0];
    cm[0][1]=a[1];
    cm[3][3]=a[2];
	
    cm[1][1]=cm[2][2]=cm[0][0];	
    cm[4][4]=cm[5][5]=cm[3][3];	
    cm[0][2]=cm[1][2]=cm[0][1];
    cm[2][0]=cm[2][1]=cm[1][0]=cm[0][1];

   } else if (ns==5) {
    /* hexagonal */
     
    if (hextype == 1) {
      /* VTI */
      fprintf(stderr, "vertical transverse isotropy\n");

      cm[2][2]=a[0];
      cm[1][2]=a[1];
      cm[0][1]=a[2];
      cm[3][3]=a[3];
      cm[5][5]=a[4];

      cm[0][0]=cm[1][1]=2.0*cm[5][5] + cm[0][1];
      cm[0][2]=cm[2][0]=cm[2][1]=cm[1][2];
      cm[1][0]=cm[0][1];
      cm[4][4]=cm[3][3];
    }

    else if (hextype == 2) {
      /* HTI */
      fprintf(stderr, "horizontal transverse isotropy\n");

      cm[0][0]=a[0];
      cm[2][2]=a[1];
      cm[0][1]=a[2];
      cm[5][5]=a[3];
      cm[3][3]=a[4];

      cm[1][2]=cm[2][1]=cm[2][2] - 2.0*cm[3][3];
      cm[0][2]=cm[1][0]=cm[2][0]=cm[0][1];
      cm[1][1]=cm[2][2];
      cm[4][4]=cm[5][5];

    }
  }

    else if (ns==6){
      /* tetragonal */
      cm[0][0]=a[0];
      cm[2][2]=a[1];
      cm[1][2]=a[2];
      cm[0][1]=a[3];
      cm[3][3]=a[4];
      cm[5][5]=a[5];

      cm[1][1]=cm[0][0];
      cm[0][2]=cm[2][0]=cm[1][2];
      cm[1][0]=cm[0][1];
      cm[2][1]=cm[1][2];
      cm[4][4]=cm[3][3];
  
  } 
  
    else if (ns==9){
      /* orthorhombic */
      cm[0][0]=a[0];
      cm[1][1]=a[1];
      cm[2][2]=a[2];
      cm[1][2]=a[3];
      cm[0][2]=a[4];
      cm[0][1]=a[5];
      cm[3][3]=a[6];
      cm[4][4]=a[7];
      cm[5][5]=a[8];

      cm[2][0]=cm[0][2];
      cm[1][0]=cm[0][1];
      cm[2][1]=cm[1][2];
  } 
    else {
    err("given elatic moduli does not fit given ns");
  }

  /* stiffness tensor calculation */
  c= (double ****) malloc(sizeof(double ***)*3);
  for (i=0; i<3; ++i)
    c[i]=alloc3double(3,3,3);
  stiffness (c,  cm);
       
  /* alloc workspace to solve for eigenvalues and eigenfunctions */
  e= (double **) malloc(8*sizeof(double *));
  for (k=0;  k<8; ++k)
    e[k] = alloc1double(irk[k]*irk[k]);
  
  gamma= (double **) malloc(8*sizeof(double *));
  for (k=0;  k<8; ++k)
    gamma[k] = alloc1double(irk[k]*irk[k]);
      
  /* filling matrix e */
  for (k=0; k<8; ++k)
    e_fill(e[k], itab, ltab, mtab, ntab, 
	   r, d1, d2, d3, rho, shape, k, irk);
    
  /* filling matrix gamma  */
  for (k=0; k<8; ++k)
    gamma_fill(gamma[k], itab, ltab, mtab, 
	       ntab, r, d1, d2, d3, c, shape, k, irk);
      
  
  fprintf(stderr, "starting eigenvalues calculation\n"); 
  /*-------------------------------------------------------------*/
  /*--------- solve the generalized eigenvalue problem ----------*/
  /*-------------------------------------------------------------*/  
  w= (double **) malloc(sizeof(double *)*8);
  itype=1;
  jobz='V';
  uplo='U';
  for (k=0; k<8; ++k){
    w[k] =alloc1double(irk[k]);
    lda=ldb=irk[k]; 
    order=irk[k];  
    lwork=MAX(1, 3*order-1);
    work=alloc1double(lwork);
    /* lapack routine */
    dsygv_(&itype, &jobz, &uplo, &order, gamma[k], &lda, e[k], &ldb, 
	   w[k], work, &lwork, &info);
    free1double(work);
  }
  /*-------------------------------------------------------------*/  
  /*-------------------------------------------------------------*/
  /*-------------------------------------------------------------*/
  /* eigen vectors */
  z=alloc2double(r,r);
  irf=0;
  for (k=0; k<8; ++k){
    for (ir1=0;ir1<irf;++ir1)
      for (ir2=irf;ir2<irf+irk[k];++ir2)
	z[ir2][ir1]=0.0;
    for (ir1=irf; ir1<irf+irk[k]; ++ir1)
      for (ir2=irf; ir2<irf+irk[k]; ++ir2)
	z[ir2][ir1]=gamma[k][(ir2-irf)*irk[k]+ir1-irf]; 
    /* change the order of the array at the same time since we go
       from fortran array to C array */
    for (ir1=irk[k]+irf;ir1<r;++ir1)
      for (ir2=irf; ir2<irf+irk[k]; ++ir2)
	z[ir2][ir1]=0.0;
    irf+=irk[k];
  }
  /* free workspace */
  for (k=0; k<8; ++k){
    free1double(gamma[k]);
    free1double(e[k]);
  }
  free(gamma);
  free(e);
  for (i=0; i<3; ++i)
    free3double(c[i]);
  free(c);
  /*fprintf(stderr, "eigenvalues calculation done\n");*/
  /* sort eigenfrequencies */
  wsort=alloc1double(r);
  wnosort=alloc1double(r);
  indice=alloc1int(r);
  for (ir1=0;ir1<r; ++ir1)
    indice[ir1]=ir1;
  for (i=0, k=0; k<8; ++k){
    for (ir1=0;ir1<irk[k];++ir1,++i){
      wsort[i]=w[k][ir1];
      wnosort[i]=w[k][ir1];
    }
  }
  dqksort(r,wsort);
  dqkisort(r,wnosort,indice);
  /* frequencies in MegaHertz */
  irf=-1;
  for (ir1=0; (ir1<r) ; ++ir1)
    if ((wsort[ir1]>0) && ((sqrt(wsort[ir1])/(2.0*PI))>0.00001))
      wsort[ir1]=sqrt(wsort[ir1])/(2.0*PI);
    else wsort[ir1]=0.0;
  
  /* output */
  dxindex(r, wsort, freqmin, &ifw);
  ++ifw;
  for (i=0, iw=ifw; i<ndata; ++i, ++iw){
    y[i]=wsort[iw];
     /* fprintf(stderr,"f%d=%i\n",i+1, (double)(y[i])); */
    fprintf(stderr,"f%d=%f\n",i+1,y[i]); 
  }

  /* write predicted frequencies to file predictedf */
  freqfile = fopen(freqs,"w");
  for (i=0, iw=ifw; i<ndata; ++i, ++iw){
    y[i]=wsort[iw];
    fprintf(freqfile,"%f\n",y[i]);
  }
  fclose(freqfile);


  /* compute dyda */
  compute_dyda(dyda,ns,hextype,r,itab,ltab,mtab,ntab,d1,d2,d3,shape,ifw,
		      ndata,z,wsort,indice);
  
  /* print out gradient of objective function */
  /* print_stuff(ns,dyda,ndata);*/
      
  /* clean workspace */
  for (i=0, k=0; k<8; ++k)
    free1double(w[k]);
  free(w);
}

/*void print_stuff(int ns, double **dyda, int ndata){*/
/*  int k,i;*/
/*  fprintf(stderr,"\nns=%d\n",ns);*/
/*  for (k=0; k<ns;k++){*/
/*    for (i=0; i<ndata; ++i){*/
/*      fprintf(stderr,"dyda[%d][%d]=%f\n",k,i,dyda[k][i]);*/
/*    }*/
/*    fprintf(stderr,"\n");*/
/*  }*/
/*}*/

void compute_dyda(double **dyda,int ns, int hextype, int r,int *itab,int *ltab,
		  int *mtab,int *ntab,double d1,double d2,double d3,
		  int shape,int ifw,int ndata,double **z,double *wsort,
		  int *indice){
  int i;
  int iw;

  double ****dc_c11;
  double ****dc_c22;
  double ****dc_c12;
  double ****dc_c13;
  double ****dc_c44;
  double ****dc_c55;
  double ****dc_c33;
  double ****dc_c23;
  double ****dc_c66;
  double **dgamma_c11;
  double **dgamma_c22;
  double **dgamma_c12;
  double **dgamma_c13;
  double **dgamma_c44;
  double **dgamma_c55;
  double **dgamma_c33;
  double **dgamma_c23;
  double **dgamma_c66;
  
  /* isotropic case */  
  if (ns==2) {
    dc_c11= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c11[i]=alloc3double(3,3,3);
    dstiff_iso_c11 (dc_c11);
  
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_iso_c44 (dc_c44);
    
    /* fill dgamma_c11 */
    dgamma_c11=alloc2double(r, r);
    dgamma_fill(dgamma_c11,itab,ltab,mtab,ntab,r,d1,d2,d3,dc_c11,shape);
    for (i=0; i<3; ++i)
      free3double(dc_c11[i]);
    free(dc_c11);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44,itab,ltab,mtab,ntab,r,d1,d2,d3,dc_c44,shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c11, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
    }
    free2double(dgamma_c11);
    free2double(dgamma_c44);
  } 

    /* cubic */
    else if (ns==3) {
      
    dc_c11= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c11[i]=alloc3double(3,3,3);
    dstiff_cub_c11 (dc_c11);
  
    dc_c12= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c12[i]=alloc3double(3,3,3);
    dstiff_cub_c12 (dc_c12);
  
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_cub_c44 (dc_c44);
    
    /* fill dgamma_c11 */
    dgamma_c11=alloc2double(r, r);
    dgamma_fill(dgamma_c11,itab,ltab,mtab,ntab,r,d1,d2,d3,dc_c11,shape);
    for (i=0; i<3; ++i)
      free3double(dc_c11[i]);
    free(dc_c11);
    
    /* fill dgamma_c12 */
    dgamma_c12=alloc2double(r, r);
    dgamma_fill(dgamma_c12, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c12, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c12[i]);
    free(dc_c12);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c44, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c11, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c12, z, indice[iw], r);
      dyda[2][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
    }
    free2double(dgamma_c11);
    free2double(dgamma_c12);
    free2double(dgamma_c44);

 
  } 

    else if (ns==5) {
      /* hexagonal */

    if (hextype == 1){
      /* VTI */
    fprintf(stderr, "hextype=%d\n", hextype);
  

    dc_c33= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c33[i]=alloc3double(3,3,3);
    dstiff_vti_c33 (dc_c33);
  
    dc_c23= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c23[i]=alloc3double(3,3,3);
    dstiff_vti_c23 (dc_c23);
  
    dc_c12= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c12[i]=alloc3double(3,3,3);
    dstiff_vti_c12 (dc_c12);
    
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_vti_c44 (dc_c44);
    
    dc_c66= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c66[i]=alloc3double(3,3,3);
    dstiff_vti_c66 (dc_c66);
    
    /* fill dgamma_c33 */
    dgamma_c33=alloc2double(r, r);
    dgamma_fill(dgamma_c33, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c33, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c33[i]);
    free(dc_c33);
    
    /* fill dgamma_c23 */
    dgamma_c23=alloc2double(r, r);
    dgamma_fill(dgamma_c23, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c23, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c23[i]);
    free(dc_c23);
    
    /* fill dgamma_c12 */
    dgamma_c12=alloc2double(r, r);
    dgamma_fill(dgamma_c12, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c12, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c12[i]);
    free(dc_c12);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c44, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* fill dgamma_c66 */
    dgamma_c66=alloc2double(r, r);
    dgamma_fill(dgamma_c66, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c66, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c66[i]);
    free(dc_c66);
    
    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c33, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c23, z, indice[iw], r);
      dyda[2][i]=dfdp(wsort[iw], dgamma_c12, z, indice[iw], r);
      dyda[3][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
      dyda[4][i]=dfdp(wsort[iw], dgamma_c66, z, indice[iw], r);
      
   
    }
    free2double(dgamma_c33);
    free2double(dgamma_c23);
    free2double(dgamma_c12);
    free2double(dgamma_c44);
    free2double(dgamma_c66);
    
      }
   

    else if (hextype == 2) {
      /* HTI */
      fprintf(stderr, "hextype=%d\n", hextype);
   
    dc_c11= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c11[i]=alloc3double(3,3,3);
    dstiff_hti_c11 (dc_c11);

    dc_c33= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c33[i]=alloc3double(3,3,3);
    dstiff_hti_c33 (dc_c33);
  
    dc_c12= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c12[i]=alloc3double(3,3,3);
    dstiff_hti_c12 (dc_c12);
    
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_hti_c44 (dc_c44);
    
    dc_c66= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c66[i]=alloc3double(3,3,3);
    dstiff_hti_c66 (dc_c66);

    /* fill dgamma_c11 */
    dgamma_c11=alloc2double(r, r);
    dgamma_fill(dgamma_c11, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c33, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c11[i]);
    free(dc_c11);

    /* fill dgamma_c33 */
    dgamma_c33=alloc2double(r, r);
    dgamma_fill(dgamma_c33, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c33, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c33[i]);
    free(dc_c33);
      
    /* fill dgamma_c12 */
    dgamma_c12=alloc2double(r, r);
    dgamma_fill(dgamma_c12, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c12, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c12[i]);
    free(dc_c12);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c44, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* fill dgamma_c66 */
    dgamma_c66=alloc2double(r, r);
    dgamma_fill(dgamma_c66, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c66, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c66[i]);
    free(dc_c66);

    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c11, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c33, z, indice[iw], r);
      dyda[2][i]=dfdp(wsort[iw], dgamma_c12, z, indice[iw], r);
      dyda[3][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
      dyda[4][i]=dfdp(wsort[iw], dgamma_c66, z, indice[iw], r);

    }
    free2double(dgamma_c11);
    free2double(dgamma_c33);
    free2double(dgamma_c12);
    free2double(dgamma_c44);
    free2double(dgamma_c66);

    }

 }

  else if (ns==6) {
    /* tetragonal */
    dc_c11= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c11[i]=alloc3double(3,3,3);
    dstiff_tetra_c11 (dc_c11);

    dc_c33= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c33[i]=alloc3double(3,3,3);
    dstiff_tetra_c33 (dc_c33);
  
    dc_c23= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c23[i]=alloc3double(3,3,3);
    dstiff_tetra_c23 (dc_c23);
  
    dc_c12= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c12[i]=alloc3double(3,3,3);
    dstiff_tetra_c12 (dc_c12);
    
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_tetra_c44 (dc_c44);
    
    dc_c66= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c66[i]=alloc3double(3,3,3);
    dstiff_tetra_c66 (dc_c66);
    /* fill dgamma_c11 */
    dgamma_c11=alloc2double(r, r);
    dgamma_fill(dgamma_c11, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c11, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c11[i]);
    free(dc_c11);
    
    /* fill dgamma_c33 */
    dgamma_c33=alloc2double(r, r);
    dgamma_fill(dgamma_c33, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c33, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c33[i]);
    free(dc_c33);
    
    /* fill dgamma_c23 */
    dgamma_c23=alloc2double(r, r);
    dgamma_fill(dgamma_c23, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c23, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c23[i]);
    free(dc_c23);
    
    /* fill dgamma_c12 */
    dgamma_c12=alloc2double(r, r);
    dgamma_fill(dgamma_c12, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c12, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c12[i]);
    free(dc_c12);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c44, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* fill dgamma_c66 */
    dgamma_c66=alloc2double(r, r);
    dgamma_fill(dgamma_c66, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c66, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c66[i]);
    free(dc_c66);
    
    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c11, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c33, z, indice[iw], r);
      dyda[2][i]=dfdp(wsort[iw], dgamma_c23, z, indice[iw], r);
      dyda[3][i]=dfdp(wsort[iw], dgamma_c12, z, indice[iw], r);
      dyda[4][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
      dyda[5][i]=dfdp(wsort[iw], dgamma_c66, z, indice[iw], r);
    }
    free2double(dgamma_c11);
    free2double(dgamma_c33);
    free2double(dgamma_c23);
    free2double(dgamma_c12);
    free2double(dgamma_c44);
    free2double(dgamma_c66);
  } 

  else if (ns==9) {
    /* orthorhombic */
    dc_c11= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c11[i]=alloc3double(3,3,3);
    dstiff_orth_c11 (dc_c11);

    dc_c22= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c22[i]=alloc3double(3,3,3);
    dstiff_orth_c22 (dc_c22);

    dc_c33= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c33[i]=alloc3double(3,3,3);
    dstiff_orth_c33 (dc_c33);
  
    dc_c23= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c23[i]=alloc3double(3,3,3);
    dstiff_orth_c23 (dc_c23);
    
    dc_c13= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c13[i]=alloc3double(3,3,3);
    dstiff_orth_c13 (dc_c13);
  
    dc_c12= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c12[i]=alloc3double(3,3,3);
    dstiff_orth_c12 (dc_c12);
    
    
    dc_c44= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c44[i]=alloc3double(3,3,3);
    dstiff_orth_c44 (dc_c44);
    
    dc_c55= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c55[i]=alloc3double(3,3,3);
    dstiff_orth_c55 (dc_c55);
    
    dc_c66= (double ****) malloc(sizeof(double ***)*3);
    for (i=0; i<3; ++i)
      dc_c66[i]=alloc3double(3,3,3);
    dstiff_orth_c66 (dc_c66);
    
    /* fill dgamma_c11 */
    dgamma_c11=alloc2double(r, r);
    dgamma_fill(dgamma_c11, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c11, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c11[i]);
    free(dc_c11);
    /* fill dgamma_c22 */
    dgamma_c22=alloc2double(r, r);
    dgamma_fill(dgamma_c22, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c22, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c22[i]);
    free(dc_c22);
    
    /* fill dgamma_c33 */
    dgamma_c33=alloc2double(r, r);
    dgamma_fill(dgamma_c33, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c33, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c33[i]);
    free(dc_c33);
    
    /* fill dgamma_c23 */
    dgamma_c23=alloc2double(r, r);
    dgamma_fill(dgamma_c23, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c23, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c23[i]);
    free(dc_c23);
    
    /* fill dgamma_c23 */
    dgamma_c13=alloc2double(r, r);
    dgamma_fill(dgamma_c13, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c13, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c13[i]);
    free(dc_c13);
    
    /* fill dgamma_c12 */
    dgamma_c12=alloc2double(r, r);
    dgamma_fill(dgamma_c12, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c12, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c12[i]);
    free(dc_c12);
    
    /* fill dgamma_c44 */
    dgamma_c44=alloc2double(r, r);
    dgamma_fill(dgamma_c44, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c44, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c44[i]);
    free(dc_c44);

    /* fill dgamma_c55 */
    dgamma_c55=alloc2double(r, r);
    dgamma_fill(dgamma_c55, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c55, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c55[i]);
    free(dc_c55);

    /* fill dgamma_c66 */
    dgamma_c66=alloc2double(r, r);
    dgamma_fill(dgamma_c66, itab, ltab, mtab, ntab, r, d1, d2, d3, dc_c66, shape);
    for (i=0; i<3; ++i)
      free3double(dc_c66[i]);
    free(dc_c66);
    
    /* gradiant of the objective function */
    for (iw=ifw,i=0; i<ndata; ++i, ++iw){
      dyda[0][i]=dfdp(wsort[iw], dgamma_c11, z, indice[iw], r);
      dyda[1][i]=dfdp(wsort[iw], dgamma_c22, z, indice[iw], r);
      dyda[2][i]=dfdp(wsort[iw], dgamma_c33, z, indice[iw], r);
      dyda[3][i]=dfdp(wsort[iw], dgamma_c23, z, indice[iw], r);
      dyda[4][i]=dfdp(wsort[iw], dgamma_c13, z, indice[iw], r);
      dyda[5][i]=dfdp(wsort[iw], dgamma_c12, z, indice[iw], r);
      dyda[6][i]=dfdp(wsort[iw], dgamma_c44, z, indice[iw], r);
      dyda[7][i]=dfdp(wsort[iw], dgamma_c55, z, indice[iw], r);
      dyda[8][i]=dfdp(wsort[iw], dgamma_c66, z, indice[iw], r);
    }
    free2double(dgamma_c11); 
    free2double(dgamma_c22);
    free2double(dgamma_c33);
    free2double(dgamma_c23);
    free2double(dgamma_c13);
    free2double(dgamma_c12);
    free2double(dgamma_c44);
    free2double(dgamma_c55);
    free2double(dgamma_c66);
  }
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
	       *doublefact(l-1)*doublefact(m-1)*doublefact(n-1)/doublefact(l+m+n+3);
    /* rp shape */
    default: return 8.0/((l+1)*(m+1)*(n+1))*pow(d1, l+1)*pow(d2,m+1)*pow(d3, n+1);
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
/* isotropic */
void dstiff_iso_c11(double ****c){
  double **cm;
  int i,j;
  
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][0]=10.0;
  cm[1][1]=cm[2][2]=cm[0][0];
  cm[0][1]=cm[0][2]=cm[1][2]=cm[0][0];
  cm[1][0]=cm[2][0]=cm[2][1]=cm[0][0];
  stiffness (c,cm);
}
void dstiff_iso_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
  cm[4][4]=cm[5][5]=cm[3][3];	
  cm[0][1]=cm[0][2]=cm[1][2]=-2.0*cm[3][3];
  cm[1][0]=cm[2][0]=cm[2][1]=-2.0*cm[3][3];
  stiffness (c,cm);
}

/* cubic */
void dstiff_cub_c11(double ****c){
  double **cm;
  int i,j;
  
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][0]=10.0;
  cm[1][1]=cm[2][2]=cm[0][0];	
  
  stiffness (c,cm);
}
void dstiff_cub_c12(double ****c){
  double **cm;
  int i,j;
  
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][1]=10.0;
  cm[0][2]=cm[1][2]=cm[0][1];
  cm[2][0]=cm[2][1]=cm[1][0]=cm[0][1];
  
  stiffness (c,cm);
}
void dstiff_cub_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
  cm[4][4]=cm[5][5]=cm[3][3];	
  stiffness (c,cm);
}

/* VTI */
void dstiff_vti_c33(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[2][2]=10.0;
  stiffness (c,cm);
}

void dstiff_vti_c23(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[1][2]=10.0;
  cm[0][2]=cm[2][0]=cm[2][1]=cm[1][2];
  stiffness (c,cm);
}

void dstiff_vti_c12(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][1]=10.0;
  cm[0][0]=cm[1][1]=cm[0][1];
  cm[1][0]=cm[0][1];
  stiffness (c,cm);
}

void dstiff_vti_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
  cm[4][4]=cm[3][3];
  stiffness (c,cm);
}

void dstiff_vti_c66(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[5][5]=1.0;
  cm[0][0]=cm[1][1]=2.0*cm[5][5];
  stiffness (c,cm);
}

/* HTI */
void dstiff_hti_c11(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][0]=10.0;
  stiffness (c,cm);
}

void dstiff_hti_c33(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[2][2]=10.0;
  cm[1][1]=cm[2][2];
  stiffness (c,cm);
}
  

void dstiff_hti_c12(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][1]=10.0;
  cm[0][2]=cm[1][0]=cm[2][0]=cm[0][1];
  stiffness (c,cm);
}

void dstiff_hti_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
  stiffness (c,cm);
}

void dstiff_hti_c66(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[5][5]=1.0;
  cm[4][4]=cm[5][5];
  stiffness (c,cm);
}




/* Tetragonal */
void dstiff_tetra_c11(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[1][1]=cm[0][0]=10.0;
  stiffness (c,cm);
}

void dstiff_tetra_c33(double ****c){ 
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[2][2]=10.0;
  stiffness (c,cm);
}

void dstiff_tetra_c23(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[1][2]=10.0;
  cm[0][2]=cm[2][0]=cm[2][1]=cm[1][2];
  stiffness (c,cm);
}
void dstiff_tetra_c12(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][1]=10.0;
  cm[1][0]=cm[0][1];
  stiffness (c,cm);
}

void dstiff_tetra_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
  cm[4][4]=cm[3][3];
  stiffness (c,cm);
}

void dstiff_tetra_c66(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[5][5]=1.0;
  stiffness (c,cm);
}


/* Orthorhombic */
void dstiff_orth_c11(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][0]=10.0;
  stiffness (c,cm);
}
void dstiff_orth_c22(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[1][1]=10.0;
  stiffness (c,cm);
}

void dstiff_orth_c33(double ****c){ 
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[2][2]=10.0;
  stiffness (c,cm);
}

void dstiff_orth_c23(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[1][2]=10.0;
  cm[2][1]=cm[1][2];
  stiffness (c,cm);
}

void dstiff_orth_c13(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][2]=10.0;
  cm[2][0]=cm[0][2];
  stiffness (c,cm);
}

void dstiff_orth_c12(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[0][1]=10.0;
  cm[1][0]=cm[0][1];
  stiffness (c,cm);
}

void dstiff_orth_c44(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[3][3]=1.0;
 stiffness (c,cm);
}
void dstiff_orth_c55(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[4][4]=1.0;
 stiffness (c,cm);
}

void dstiff_orth_c66(double ****c){
  double **cm;
  int i,j;
  cm=alloc2double(6,6);
  for (i=0; i<6; ++i)
    for (j=0; j<6; ++j)
      cm[i][j]=0.0;
  cm[5][5]=1.0;
  stiffness (c,cm);
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

void dgamma_fill(double **dgamma, int *itab, int *ltab, int *mtab, int *ntab, 
	   int r, double d1, double  d2, double d3, double ****dc, int shape)
{
  int ir1, ir2, i1, i2, l1, l2, m1, m2, n1, n2;
  int j1, j2;
  int l, m, n;


  
  for (ir1=0; ir1<r; ++ir1){
    for (ir2=0; ir2<r; ++ir2){
      i1=itab[ir1];
      i2=itab[ir2];
      l1=ltab[ir1];
      l2=ltab[ir2];
      m1=mtab[ir1];
      m2=mtab[ir2];
      n1=ntab[ir1];
      n2=ntab[ir2];
      dgamma[ir1][ir2]=0.0;
      if (l1>0) {
	j1=0;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-2;
	  m=m1+m2;
	  n=n1+n2;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*l1*l2*volintegral(d1, d2, d3, l, m, n, shape);
	  
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2-1;
	  m=m1+m2-1;
	  n=n1+n2;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*l1*m2*volintegral(d1, d2, d3, l, m, n, shape);
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2-1;
	  m=m1+m2;
	  n=n1+n2-1;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*l1*n2*volintegral(d1, d2, d3, l, m, n, shape);
	  
	}
      }
      if (m1>0) {
	j1=1;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-1;
	  m=m1+m2-1;
	  n=n1+n2;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*m1*l2*volintegral(d1, d2, d3, l, m, n, shape);
	  
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2;
	  m=m1+m2-2;
	  n=n1+n2; 
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*m1*m2*volintegral(d1, d2, d3, l, m, n, shape);
	  
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2;
	  m=m1+m2-1;
	  n=n1+n2-1;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*m1*n2*volintegral(d1, d2, d3, l, m, n, shape);
	
	}
      }
      if (n1>0) {
	j1=2;
	if (l2>0) {
	 j2=0;
	  l=l1+l2-1;
	  m=m1+m2;
	  n=n1+n2-1;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*n1*l2*volintegral(d1, d2, d3, l, m, n, shape);
	
	}
	if (m2>0) {
	  j2=1;
	  l=l1+l2;
	  m=m1+m2-1;
	  n=n1+n2-1;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*n1*m2*volintegral(d1, d2, d3, l, m, n, shape);
	
	}
	if (n2>0) {
	  j2=2;
	  l=l1+l2;
	  m=m1+m2;
	  n=n1+n2-2;
	  dgamma[ir1][ir2] +=dc[i1][j1][i2][j2]*n1*n2*volintegral(d1, d2, d3, l, m, n, shape);
	
	}
      }
    }
  }
}

double  dfdp(double f, double **dgammadp, double **z, int ie, int n)
{
  int i;
  double *p; 
  p=alloc1double(n);
  
  for (i=0; i<n; ++i)
    p[i]= scalarproduct(dgammadp[i], z[ie], n); /* we use dgammadp's 
						   symmetry here*/
 
  return scalarproduct(z[ie], p, n)/(8.0*PI*PI*f); 
  free1double(p);
}

double  scalarproduct(double *a, double *b, int n)
{
  int i; 
  double sp; 
  sp=0.0; 
  for (i=0; i<n; ++i)
    sp+=a[i]*b[i];
  return sp;
}


void dxindex (int nx, double ax[], double x, int *index)
/*****************************************************************************
determine index of x with respect to an array of x values
******************************************************************************
Input:
nx		number of x values in array ax
ax		array[nx] of monotonically increasing or decreasing x values
x		the value for which index is to be determined
index		index determined previously (used to begin search)

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x
******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
{
	int lower,upper,middle,step;

	/* initialize lower and upper indices and step */
	lower = *index;
	if (lower<0) lower = 0;
	if (lower>=nx) lower = nx-1;
	upper = lower+1;
	step = 1;

	/* if x values increasing */
	if (ax[nx-1]>ax[0]) {

		/* find indices such that ax[lower] <= x < ax[upper] */
		while (lower>0 && ax[lower]>x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]<=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x>=ax[middle])
				lower = middle;
			else
				upper = middle;
		}

	/* else, if not increasing */
	} else {

		/* find indices such that ax[lower] >= x > ax[upper] */
		while (lower>0 && ax[lower]<x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]>=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x<=ax[middle])
				lower = middle;
			else
				upper = middle;
		}
	}

	/* return lower index */
	*index = lower;
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

void
dqkipart (double a[], int i[], int p, int q, int *j, int *k)
/*****************************************************************************
quicksort partition (FOR INTERNAL USE ONLY):
take the value x of a random element from the subarray a[p:q] of
a[0:n-1] and rearrange indices in the subarray i[p:q] in such a way
that there exist integers j and k with the following properties:
  p <= j < k <= q, provided that p < q
  a[i[l]] <= x,  for p <= l <= j
  a[i[l]] == x,  for j < l < k
  a[i[l]] >= x,  for k <= l <= q
note that this effectively partitions the subarray with bounds
[p:q] into lower and upper subarrays with bounds [p:j] and [k:q]
******************************************************************************
Input:
a		array[p:q]
i		array[p:q] of indices to be rearranged
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
i		array[p:q] of indices rearranged
j		upper bound of lower output subarray
k		lower bound of upper output subarray
******************************************************************************
Notes:
This function is adapted from procedure partition by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int pivot,left,right,temp;
	double apivot;
	static long int seed=0L;
 
	/* choose random pivot element between p and q, inclusive */
	seed = (seed*FA+FC)%FM;
	pivot = p+(q-p)*(double)seed/(double)FM;
	if (pivot<p) pivot = p;
	if (pivot>q) pivot = q;
	apivot = a[i[pivot]];

	/* initialize left and right pointers and loop until break */
	for (left=p,right=q;;) {
		/*
		 * increment left pointer until either
		 * (1) an element greater than the pivot element is found, or
		 * (2) the upper bound of the input subarray is reached
		 */
		while (a[i[left]]<=apivot && left<q) left++;
 
		/*
		 * decrement right pointer until either
		 * (1) an element less than the pivot element is found, or
		 * (2) the lower bound of the input subarray is reached
		 */
		while (a[i[right]]>=apivot && right>p) right--;
 
		/* if left pointer is still to the left of right pointer */
		if (left<right) {

			/* exchange left and right indices */
			temp = i[left];
			i[left++] = i[right];
			i[right--] = temp;
		} 
		/* else, if pointers are equal or have crossed, break */
		else break;
	}
	/* if left pointer has not crossed pivot */
	if (left<pivot) {

		/* exchange indices at left and pivot */
		temp = i[left];
		i[left++] = i[pivot];
		i[pivot] = temp;
	}
	/* else, if right pointer has not crossed pivot */
	else if (pivot<right) {

		/* exchange indices at pivot and right */
		temp = i[right];
		i[right--] = i[pivot];
		i[pivot] = temp;
	}
	/* left and right pointers have now crossed; set output bounds */
	*j = right;
	*k = left;
}

void
dqkiinss (double a[], int i[], int p, int q)
/*****************************************************************************
quicksort insertion sort (FOR INTERNAL USE ONLY):
Sort a subarray of indices bounded by p and q so that
a[i[p]] <= a[i[p+1]] <= ... <= a[i[q]]
******************************************************************************
Input:
a		subarray[p:q] containing elements
i		subarray[p:q] containing indices to be sorted
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
i		subarray[p:q] of indices sorted
******************************************************************************
Notes:
Adapted from Sedgewick, R., 1983, Algorithms, Addison Wesley, p. 96.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int j,k,ij;
	double aij;

	for (j=p+1; j<=q; j++) {
		for (ij=i[j],aij=a[ij],k=j; k>p && a[i[k-1]]>aij; k--)
			i[k] = i[k-1];
		i[k] = ij;
	}
}

void
dqkisort (int n, double a[], int i[])
/*****************************************************************************
Sort an array of indices i[] so that 
a[i[0]] <= a[i[1]] <= ... <= a[i[n-1]]
******************************************************************************
Input:
n		number of elements in array a
a		array[n] elements
i		array[n] indices to be sorted

Output:
i		array[n] indices sorted
******************************************************************************
Notes:
n must be less than 2^NSTACK, where NSTACK is defined above.

This function is adapted from procedure quicksort by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321;
the main difference is that recursion is accomplished
explicitly via a stack array for efficiency; also, a simple
insertion sort is used to sort subarrays too small to be
partitioned efficiently.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
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
			dqkipart(a,i,p,q,&j,&k);

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
		dqkiinss(a,i,p,q);
	}
}




/* ----------------------------------------------------------- */
/* ----------------------------------------------------------- */
/* ----------------  Optimization routines  ------------------ */
/* ----------------------------------------------------------- */
/* ----------------------------------------------------------- */


/* ------------------------------------------------------------------- */
/* ---------------------------- MRQMIN ------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------- Description from: ----------------------- */
/* -------- Numerical Recipes: the art of scientific computing ------- */
/* --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling -- */
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Subroutine MRQMIN(X, Y, SIG, NDATA, A, MA, LISTA, MFIT,
 * COVAR, ALPHA, NCA, CHISQ, FUNCS, ALAMDA)
 * Levenberg-Marquardt method, attemping to reduce the value of chisq of a
 * fit between a set of NDATA points X(I), Y(I) with individual standard 
 * deviations SIG(I), and a nonlinear function dependent of MA coefficients
 * A. The array LISTA numbers the parameters A such that the first MFIT 
 * correspond to values actually being adjusted; the remaining MA-MFIT
 * parameters are held fixed at their input value. The program returns
 * current best-fit values for the MA fit parameters A and CHISQ. The 
 * arrays COVAR(NCA, NCA) and ALPHA(NCA, NCA) with physical dimension NCA
 * are used as working space during most iterations. Supply a subroutine
 * FUNCS(X, A, YFIT, DYDA, MA) that evaluates the fitting function YFIT,
 * and its derivative DYDA with respect to the fittng parameters A at X.
 * On the first call provide an initial guess for the parameters A, and 
 * set ALAMDA<0 for initialization (which then sets ALAMDA=0.001). If the
 * step succedds CHISQ becomes smaller and ALAMDA decreases by a factor 
 * of 10. If the step fails ALAMDA grows by a factor of 10. You must call
 * this routine repeatedly until convergence is achieved. Then, make one
 * final call with ALAMDA=0, so that COVAR(I,J) returns the covariance
 * matrix, and ALPHA(I, J) the curvature matrix. */
/* ------------------------------------------------------------------- */ 

/* This is the optimisation routine which is called by MAIN once every iteration */
/* mrqmin(d,r,itab,ltab,mtab,ntab,
	   irk,d1,d2,d3,rho,shape,freqmin,
	   y,sig,ndata,guess,ia,ns,covar,alpha,&chisq, hextype, formod,&alamda) */

/* ------ */
/* Inputs */
/* ------ */

/* d  = order of polynomial fit eg) 8, 10 or 12*/
/* r = dimension of problem. Related to d by  r= 3*(d+1)*(d+2)*(d+3)/6 */

/* *itab, *ltab, *mtab, *ntab = 1D array's of type int. Related to r. _tab=alloc1int(r) */
/* *irk = 1D array on int. irk = alloc1int(8)*/
/* --------------------------------------- */
/* the * in front of the variable is a pointer */
/* char *p declares p to be a pointer to char */
/* *p = 0 means "put a zero in the byte which p points to" */
/* --------------------------------------- */

/* d1, d2, d3 = dimensions of sample, diameter, height */
/* rho = density in g/cm^3 */
/* shape = cylinder, sphere or rectangular parallelepiped */
/* freqmin = minimum frequency, from param_data */

/* y = 1D array of measured frequency data*/
/* sig = 1D array of weightings - should be individual standard deviations for each freq? Is this an issue? */
/* ndata = number of data points - frequencies used from freq_data */
/* a = initial guess of parameter (cij) values. From param_data. Same size as ns */

/* ia = 1D array of length ndata with entries = 1. This numbers the parameters such that the first MFIT parameters are adjusted and the remaining parameters are held constant */
/* */

/* ma = number of coefficients. Is equivalent to ns */
/* covar = covariance matrix. Of size ns by ns (number of cijs). Initialised as identity matrix */
/* alpha = curvature matrix. Of size ns by ns (number of cijs). Initialised as identity matrix */
/* chisq = the difference between measured and predicted frequencies. Is not a "traditional" chisq */
/* hextype = differentiates between VTI and HTI symmetry in the hexagonal case */
/* (*funcs) = forward model. Calculates the frequencies based on cij values */
/* alamda = parameter from conjugate-gradient method. Starts as <0 to initialise the routine and is changed in subsequent iterations */

void mrqmin(int d,int r,int *itab,int *ltab,int *mtab,int *ntab, 
	    int *irk,double d1,double d2,double d3, 
	    double rho,int shape,float freqmin,
	    double y[],double sig[],int ndata,double a[],int ia[],
	    int ma,double **covar,double **alpha,double *chisq, int hextype, 
	    void (*funcs)(int,int,int *,int *,int *,int *, 
			   int *,double,double,double, 
			   double,int,float,int,
			  double *, double *, double **, int, int)
	    ,double *alamda)
{

  int j,k,l; /* counting measures for FOR loops */
  static int mfit; /* number of cijs that are adjusted. The remaining (ma = ns) - mfit cij values are left unchanged */
  static double ochisq,*atry,*beta,*da,**oneda;
  
  
  /* IF loop is called if almada <0. This initializes the routine and almada = -1.0 is set in main before calling MRQMIN */
  if (*alamda <0.0) {
    atry=alloc1double(ma); /* create 1D arrays the size of ma = ns = number of cijs*/
    beta=alloc1double(ma);
    da=alloc1double(ma);
    /* are arrays of size ma (number of cijs) with entries = 0 */
    
    for (mfit=0, j=0;j<ma;++j) /* set mfit = 0. Times looped through = number of independent cijs*/
      if (ia[j]) mfit++; /* if (ia[j]) is the same as if (ia[j] > 0) */
    /* WHAT GOES ON HERE? Why is mfit=4 sometimes and mfit=5 sometimes? Why does it change? What changes ia? */
    
    
    /* mfit is changed to ns regardless of the number of frequencies used (at least for isotropic sample) */
    /* shouldn't mfit only be 1 for the first freq - so only one cij is changed */
    
    /* for hexagonal symmetry mfit starts out being four and is increased to five when more freqs are added */
    /* this is what is consistent with what we see in the inverse process */
    
    
    
    oneda=alloc2double(1,mfit); /* allocate a 2-d array of doubles - row vector the size of mfit? */
    
    /* set alamda to a small positive value (0.001) after the routine has been initialized by the negative value */
    *alamda=0.001;
    
    /* compute "chisq" - need to update to formal chisq*/
    mrqcof(d,r,itab,ltab,mtab,ntab,
	   irk,d1,d2,d3, 
	   rho,shape,freqmin,y,sig,ndata,a,ia,ma,alpha,beta,chisq,hextype, funcs);
    ochisq=(*chisq); /* update chisq value */
    for (j=0;j<ma;++j) atry[j]=a[j]; /* set atry[j] equal to the inital guess of the cij value a[j]*/
  }
  /* end of IF loop - for initialization of optimization routine */
  
  
  for (j=0;j<mfit;++j) {/* mfit started = 0 and then increased in the prior IF loop, but mfit <= ns */
    for (k=0;k<mfit;++k) covar[j][k]=alpha[j][k]; /* set components in covariance matrix as equal to curvature matrix. covar and alpha started as identity matrices. Must be changed in mrqcof otherwise this line woudln't do anything*/
    covar[j][j]=alpha[j][j]*(1.0+(*alamda)); /* change the main diagonals */
    oneda[j][0]=beta[j]; /* update oneda values by beta - which is changed (?) by mrqcof*/
  }
  
  gaussj(covar,mfit,oneda,1); 
  for (j=0; j<mfit; ++j) da[j]=oneda[j][0];
  
  if (*alamda==0.0) { /* this is the stopping criteria - but how do we get alamda = 0?*/
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    free2double(oneda);
    free1double(da);
    free1double(beta);
    free1double(atry);
    return;
  }
  
  for (j=-1,l=0;l<ma;++l)
    if (ia[l]) atry[l]=a[l]+da[++j];
  /* compute "chisq" - need to update to formal chisq*/
  mrqcof(d,r,itab,ltab,mtab,ntab, 
	   irk,d1,d2,d3, 
	 rho,shape,freqmin,y,sig,ndata,atry,ia,ma,covar,da,chisq, hextype, funcs);
  
  if (*chisq < ochisq) {/* if step succeeds value of chisq decreases: ochisq < chisq */
    *alamda *= 0.1; /* decrease alamda by a factor of ten */
    ochisq=(*chisq); /* update chisq with the new, reduced, value */
    for (j=0;j<mfit;++j){
      for (k=0; k<mfit;++k) alpha[j][k]=covar[j][k];
      beta[j]=da[j];
    }
    for (l=0;l<ma;++l) a[l]=atry[l];
  } else {/* else step does not succeed and chisq increases*/
    *alamda *=10.0; /* increase alamda by a factor of ten, and then loop back and try again*/
    *chisq=ochisq;
    
  }
}


# -------------------------------------------------------------------
# ---------------------------- COVSRT -------------------------------
# -------------------------------------------------------------------
# ------------------------- Description from: -----------------------
# -------- Numerical Recipes: the art of scientific computing -------
# --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling --
# -------------------------------------------------------------------

# subroutine COVSRT(COVAR, NCVM, MA, LISTA, MFIT)
# Given the covariance matrix COVAR of a fir for MFIT of MA total 
# parameters, and their ordering LISTA(I), repack the covariance matrix 
# to the true order of the parameters. Elements associated with fixed
# parameters will be zero. NCVM is the physical dimension of COVAR

def covsrt(covar, ma, ia, mfit):
	for i in range(mfit, ma):
		for j in range(i):
			covar[i][j] = covar[j][i] = 0.0
	k = mfit - 1;
	for j in range(ma-1, -1, -1):
		if ia[j]:
			for i in range(ma):
				covar[i][k], covar[i][j] = covar[i][j], covar[i][k]
			for i in range(ma):
				covar[k][i], covar[j][i] = covar[j][i], covar[k][i]
		k -= 1;



/* ------------------------------------------------------------------- */
/* ---------------------------- MRQCOF ------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------- Description from: ----------------------- */
/* -------- Numerical Recipes: the art of scientific computing ------- */
/* --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling -- */
/* ------------------------------------------------------------------- */

/* Subroutine MRQCOF(X, Y, SIG, NDATA, A, MA, LISTA, MFIT, ALPHA, BETA,
 * NALP, CHISQ, FUNCS)
 * Used by MRQMIN to evaluate the linearized fitting matrix ALPHA, and 
 * vector BETA from (14.4.8) */ 

/* called by mrqmin each iteration. Computes "chisq" */
void mrqcof(int d,int r,int *itab,int *ltab,int *mtab,int *ntab, 
	    int *irk,double d1,double d2,double d3, 
	    double rho,int shape,float freqmin,
	    double y[],double sig[],int ndata,double a[],
	    int ia[],int ma,double **alpha,double beta[],double *chisq, int hextype,
	    void (*funcs)(int,int,int *,int *,int *,int *, 
			   int *,double,double,double, 
			   double,int,float,int,
			  double *, double *, double **, int, int))
{
  int i,j,k,l,m,mfit=0;
  double *ymod,wt,sig2i,dy,**dyda;
  
  dyda=alloc2double(ndata,ma); /* creates a 2D array of type double and size ndata and ma*/
  ymod=alloc1double(ndata);
  for (j=0;j<ma;++j)
    if (ia[j]) mfit++;
  for (j=0;j<mfit;++j) {
    for (k=0;k<j;++k) alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  *chisq=0.0;
  (*funcs)(d,r,itab,ltab,mtab,ntab, 
	   irk,d1,d2,d3, 
	   rho,shape,freqmin,ndata,a,ymod,dyda,ma, hextype); 
  for (i=0; i<ndata; ++i) {
    sig2i=(sig[i]*sig[i]);
    dy=y[i]-ymod[i];
    for (j=-1,l=0; l<ma; ++l) {
      if (ia[l]) {
	wt=dyda[l][i]*sig2i;
	for (j++,k=-1,m=0;m<l;m++)
	  if (ia[m]) alpha[j][++k] += wt*dyda[m][i];
	beta[j] +=dy*wt;
      }
    }
    *chisq +=dy*dy*sig2i;
  }
  
  /* chisq prints from here */
  fprintf(stderr,"chisq=%f\n\n",100.0*(*chisq));
  for (j=1;j<mfit; ++j)
    for (k=0; k<(j-1); ++k) alpha[k][j]=alpha[j][k];
  free2double(dyda);
  free1double(ymod);
}


/* ------------------------------------------------------------------- */
/* ---------------------------- GAUSSJ ------------------------------- */
/* ------------------------------------------------------------------- */
/* ------------------------- Description from: ----------------------- */
/* -------- Numerical Recipes: the art of scientific computing ------- */
/* --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling -- */
/* ------------------------------------------------------------------- */

/* Subroutine GAUSSJ(A, N, NP, B, M, MP)
 * Linear equation solution by Gauss-Jordan elimination, equation (2.1.1)
 * above. A is an input matrix of N by N elements, stored in an array of
 * physical dimensions NP by NP. B is an input matrix of N by N containing
 * the M right-hand side vectors, stored in an array of physical 
 * dimensions Np by MP. On output, A is replaced by its matrix inverse,
 * and B is replaced by the corresponding set of solution vectors. */

void gaussj(double **a,int n,double **b,int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll,*ivector();
	double big,dum,pivinv,swap;
	void nrerror();

	indxc=alloc1int(n);
	indxr=alloc1int(n);
	ipiv=alloc1int(n);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++):
				a[irow][l], a[icol][l] = a[icol][l], a[irow][l]
			for (l=0;l<m;l++):
				b[irow][l], b[icol][l] = b[icol][l], b[irow][l]
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++):
				a[k][indxr[l]], a[k][indxc[l]] = a[k][indxc[l]], a[k][indxr[l]]
	}
	free1int(ipiv);
	free1int(indxr);
	free1int(indxc);
}






void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void _exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	_exit(1);
}
