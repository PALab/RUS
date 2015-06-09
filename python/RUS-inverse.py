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

import sys
import numpy
import rus
import rus_parser as p
from math import fabs

NSTACK = 50     # maximum sort length is 2^NSTACK
NSMALL = 7      # size of array for which insertion sort is fast
FM = 7875       # constants used to generate random pivots
FA = 211
FC = 1663

# send to parser
args = p.inverse_parser(sys.argv)

# print out initial guess
for is_1 in args.a:
    print(is_1)
    is_1 = is_1 / 100

d  = args.order
d1 = args.d1 / 2.0  # half sample dimensions are used in calculations
d2 = args.d2 / 2.0
d3 = args.d3 / 2.0 

# dimension of the problem
r = 3 *(d + 1) * (d + 2) * (d + 3) / 6;

# matrices of the eigenvalue problem and for the function gradiant
measurement = "freq_data"  # CHANGE THIS LINE TO APPROPRIATE DIRECTORY

# get measured frequencies from file
f = open(measurement, "rU")
nfreq = int(f.readline())
print('nfreq = ' + str(nfreq));
freq   = []
weight = []
for i in range(nfreq):
    line = f.readline()
    nums = line.split(None, 1)
    if len(nums) != 2:
        print('Could not parse some lines in the frequency file: ' + measurement)
        f.close()
        sys.exit(-1)
    freq.append(float(nums[0]))
    weight.append(float(nums[1]))
    print('freq = ' + str(freq[-1]));
f.close()

if len(freq) != nfreq:
    print('Unexpected number of frequencies.')
    print('Expected ' + str(nfreq) + ' but read ' + str(len(freq)))
    sys.exit(-1)

# relationship between ir and l,m,n - filling tables
itab, ltab, mtab, ntab, irk = rus.index_relationship(d, r)

ifr1 = rus.xindex(freq, args.freqmin, 0)
ifr2 = rus.xindex(freq, args.freqmax, nfreq)
# add one to ifr2 so the subtraction ifr2 - ifr1 gives the number of frequencies ndata
ifr2 += 1
# ndata = number of frequencies used for inversion - first line from freq_data
ndata = ifr2 - ifr1
print('ndata=' + str(ndata))

y   = [0.0 for i in range(ndata)]
sig = [0.0 for i in range(ndata)]
ia  = [1   for i in range(ndata)]

i = 0
for ifr in range(ifr1, ifr2):
    # set y equal to freq_data values. Therefore, y is an array of length ndata
    y[i] = freq[ifr]
    # set sig equal to the weightings from the second column of freq_data. Is an array of length ndata
    sig[i] = weight[ifr]
    i += 1
    # sig is not a true 'error' but the weighting from freq_data
    # is multipled by the difference in the misfit function/chisq so modes with a zero weighting are ignored

covar  = numpy.identity(args.ns)
alpha  = numpy.identity(args.ns)
chisq  =  0.0
alamda = -1.0
niter  =  100 # set number of iterations
for i in range(niter):
    mrqmin(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,y,sig,ndata,args.a,ia,ns,covar,alpha,&chisq,hextype,formod,alamda)
    print('iter #' + str(iter))
    for is_1 in range(args.ns): # ns = dimension of symmetry
        print(str(100 * guess[is_1])) # print estimated cij values

print()
print('This calculation can be executed again with the following command:')
print('python {} --order {} --shape {} --ns {} --hextype {} --d1 {} --d2 {} --d3 {} --rho {} --freqmin {} --freqmax {} {}'.format(sys.argv[0], args.order, args.shape, args.ns, args.hextype, args.d1, args.d2, args.d3, args.rho, args.freqmin, args.freqmax, ' '.join(('--' + s + ' ' + str(v)) for k,v in args.a.iteritems())))



def formod(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,ndata,a,y,dyda,ns,hextype):

    freqs = 'predictedf' # CHANGE THIS LINE TO APPROPRIATE DIRECTORY

    cm = rus.make_cm(a,hextype)
    e = rus.e_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, args.rho, args.shape, irk)
    gamma = rus.gamma_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, cm, args.shape, irk)
      
  
    print('starting eigenvalues calculation')
    #/*-------------------------------------------------------------*/
    #/*--------- solve the generalized eigenvalue problem ----------*/
    #/*-------------------------------------------------------------*/  
    w = []
    for k in range(8):
        # lapack routine
        a, w_temp, info = lapack.dsygv(gamma[k], e[k], itype=1, jobz='V', uplo='U');  
        w.append(w_temp)

    wsort = scipy.zeros(r)
    i = 0
    for k in range(8):
        for ir1 in range(irk[k]):
            wsort[i] = w[k][ir1]
            i += 1
    wsort.sort()

    i = 0
    ir1 = 0
    while ir1 < args.nfreq:
        if ((wsort[i]>0) and ((sqrt(wsort[i])/(2.0*scipy.pi))>0.00001)):
            ir1 += 1
            print(" f%d = %f" % (ir1, 1000000*sqrt(wsort[i])/(2.0*scipy.pi)))
        i += 1

    #/*-------------------------------------------------------------*/  
    #/*-------------------------------------------------------------*/
    #/*-------------------------------------------------------------*/
    #/* eigen vectors */
    z=alloc2double(r,r);
  irf=0;
  for (k=0; k<8; ++k){
    for (ir1=0;ir1<irf;++ir1)
      for (ir2=irf;ir2<irf+irk[k];++ir2)
	z[ir2][ir1]=0.0;
    for (ir1=irf; ir1<irf+irk[k]; ++ir1)
      for (ir2=irf; ir2<irf+irk[k]; ++ir2)
	z[ir2][ir1]=gamma[k][(ir2-irf)*irk[k]+ir1-irf]; 
    #/* change the order of the array at the same time since we go
    #   from fortran array to C array */
    for (ir1=irk[k]+irf;ir1<r;++ir1)
      for (ir2=irf; ir2<irf+irk[k]; ++ir2)
	z[ir2][ir1]=0.0;
    irf+=irk[k];
  }
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
  ifw = rus.xindex(wsort, freqmin, ifw);
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
    p[i]= scalarproduct(dgammadp[i], z[ie], n); #/* we use dgammadp's 
						   #symmetry here*/
 
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


void dqkpart (double a[], int p, int q, int *j, int *k)
{
	int pivot,left,right;
	double apivot,temp;
	static long int seed=0L;
 
	#/* choose random pivot element between p and q, inclusive */
	seed = (seed*FA+FC)%FM;
	pivot = p+(q-p)*(double)seed/(double)FM;
	if (pivot<p) pivot = p;
	if (pivot>q) pivot = q;
	apivot = a[pivot];

	#/* initialize left and right pointers and loop until break */
	for (left=p,right=q;;) {
		#/*
		# * increment left pointer until either
		# * (1) an element greater than the pivot element is found, or
		# * (2) the upper bound of the input subarray is reached
		# */
		while (a[left]<=apivot && left<q) left++;
 
		#/*
		# * decrement right pointer until either
		# * (1) an element less than the pivot element is found, or
		# * (2) the lower bound of the input subarray is reached
		# */
		while (a[right]>=apivot && right>p) right--;
 
		#/* if left pointer is still to the left of right pointer */
		if (left<right) {
			#/* exchange left and right elements */
			temp = a[left];
			a[left++] = a[right];
			a[right--] = temp;
		} 
		#/* else, if pointers are equal or have crossed, break */
		else break;
	}
	#/* if left pointer has not crossed pivot */
	if (left<pivot) {

		#/* exchange elements at left and pivot */
		temp = a[left];
		a[left++] = a[pivot];
		a[pivot] = temp;
	}
	#/* else, if right pointer has not crossed pivot */
	else if (pivot<right) {

		#/* exchange elements at pivot and right */
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



#/* ----------------------------------------------------------- */
#/* ----------------------------------------------------------- */
#/* ----------------  Optimization routines  ------------------ */
#/* ----------------------------------------------------------- */
#/* ----------------------------------------------------------- */
#
#
#/* ------------------------------------------------------------------- */
#/* ---------------------------- MRQMIN ------------------------------- */
#/* ------------------------------------------------------------------- */
#/* ------------------------- Description from: ----------------------- */
#/* -------- Numerical Recipes: the art of scientific computing ------- */
#/* --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling -- */
#/* ------------------------------------------------------------------- */
#
#/* ------------------------------------------------------------------- */
#/* Subroutine MRQMIN(X, Y, SIG, NDATA, A, MA, LISTA, MFIT,
# * COVAR, ALPHA, NCA, CHISQ, FUNCS, ALAMDA)
# * Levenberg-Marquardt method, attemping to reduce the value of chisq of a
# * fit between a set of NDATA points X(I), Y(I) with individual standard 
# * deviations SIG(I), and a nonlinear function dependent of MA coefficients
# * A. The array LISTA numbers the parameters A such that the first MFIT 
# * correspond to values actually being adjusted; the remaining MA-MFIT
# * parameters are held fixed at their input value. The program returns
# * current best-fit values for the MA fit parameters A and CHISQ. The 
# * arrays COVAR(NCA, NCA) and ALPHA(NCA, NCA) with physical dimension NCA
# * are used as working space during most iterations. Supply a subroutine
# * FUNCS(X, A, YFIT, DYDA, MA) that evaluates the fitting function YFIT,
# * and its derivative DYDA with respect to the fittng parameters A at X.
# * On the first call provide an initial guess for the parameters A, and 
# * set ALAMDA<0 for initialization (which then sets ALAMDA=0.001). If the
# * step succedds CHISQ becomes smaller and ALAMDA decreases by a factor 
# * of 10. If the step fails ALAMDA grows by a factor of 10. You must call
# * this routine repeatedly until convergence is achieved. Then, make one
# * final call with ALAMDA=0, so that COVAR(I,J) returns the covariance
# * matrix, and ALPHA(I, J) the curvature matrix. */
#/* ------------------------------------------------------------------- */ 
#
#/* This is the optimisation routine which is called by MAIN once every iteration */
#/* mrqmin(d,r,itab,ltab,mtab,ntab,
#	   irk,d1,d2,d3,rho,shape,freqmin,
#	   y,sig,ndata,guess,ia,ns,covar,alpha,&chisq, hextype, formod,&alamda) */
#
#/* ------ */
#/* Inputs */
#/* ------ */
#
#/* d  = order of polynomial fit eg) 8, 10 or 12*/
#/* r = dimension of problem. Related to d by  r= 3*(d+1)*(d+2)*(d+3)/6 */
#
#/* *itab, *ltab, *mtab, *ntab = 1D array's of type int. Related to r. _tab=alloc1int(r) */
#/* *irk = 1D array on int. irk = alloc1int(8)*/
#/* --------------------------------------- */
#/* the * in front of the variable is a pointer */
#/* char *p declares p to be a pointer to char */
#/* *p = 0 means "put a zero in the byte which p points to" */
#/* --------------------------------------- */
#
#/* d1, d2, d3 = dimensions of sample, diameter, height */
#/* rho = density in g/cm^3 */
#/* shape = cylinder, sphere or rectangular parallelepiped */
#/* freqmin = minimum frequency, from param_data */
#
#/* y = 1D array of measured frequency data*/
#/* sig = 1D array of weightings - should be individual standard deviations for each freq? Is this an issue? */
#/* ndata = number of data points - frequencies used from freq_data */
#/* a = initial guess of parameter (cij) values. From param_data. Same size as ns */
#
#/* ia = 1D array of length ndata with entries = 1. This numbers the parameters such that the first MFIT parameters are adjusted and the remaining parameters are held constant */
#/* */
#
#/* ma = number of coefficients. Is equivalent to ns */
#/* covar = covariance matrix. Of size ns by ns (number of cijs). Initialised as identity matrix */
#/* alpha = curvature matrix. Of size ns by ns (number of cijs). Initialised as identity matrix */
#/* chisq = the difference between measured and predicted frequencies. Is not a "traditional" chisq */
#/* hextype = differentiates between VTI and HTI symmetry in the hexagonal case */
#/* (*funcs) = forward model. Calculates the frequencies based on cij values */
#/* alamda = parameter from conjugate-gradient method. Starts as <0 to initialise the routine and is changed in subsequent iterations */
#
def mrqmin(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,y,sig,ndata,a,ia,ma,covar,alpha,chisq,hextype,funcs,alamda):
    # Loop is called if almada <0.
    # This initializes the routine.
    # Sets almada = -1.0 in main before calling MRQMIN.
    if alamda < 0.0:
        # create 1D arrays the size of ma = ns = number of cijs
        atry = numpy.zeros(ma)
        beta = numpy.zeros(ma)
        da   = numpy.zeros(ma)

        # mfit is the number of cijs that are adjusted.
        # The remaining (ma = ns) - mfit cij values are left unchanged.
        # Set mfit = 0. Times looped through = number of independent cijs.
        mfit = 0
        for j in range(ma):
            if ia[j] > 0:
                mfit += 1

        # allocate a 2-d array of doubles - row vector the size of mfit?
        oneda = numpy.zeros((1,mfit))

        # Set alamda to a small positive value (0.001)
        # after the routine has been initialized by the negative value
        alamda = 0.001

        # Compute "chisq" - need to update to formal chisq.
        mrqcof(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,y,sig,ndata,a,ia,ma,alpha,beta,chisq,hextype,funcs)

        # update chisq value
        ochisq = chisq
        for j in range(ma):
            # set atry[j] equal to the inital guess of the cij value a[j]
            atry[j] = a[j]
  
    # mfit started = 0 and then increased in the prior if loop, but mfit <= ns
    for j in range(mfit):
        # Set components in covariance matrix as equal to curvature matrix.
        # covar and alpha started as identity matrices.
        # Must be changed in mrqcof() otherwise this line would not do anything.
        for k in range(mfit):
            covar[j][k] = alpha[j][k]
        # change the main diagonals
        covar[j][j] = alpha[j][j] * (1.0 + alamda)
        # update oneda values by beta - which is changed (?) by mrqcof()
        oneda[j][0] = beta[j]
  
    gaussj(covar,mfit,oneda,1)
    for j in range(mfit):
        da[j] = oneda[j][0]

    # This is the stopping criteria - but how do we get alamda = 0.0?
    if alamda == 0.0:
        covsrt(covar,ma,ia,mfit)
        covsrt(alpha,ma,ia,mfit)
        return

    j = -1
    for l in range(ma):
        if ia[l] != 0:
            j += 1
            atry[l] = a[l] + da[j]

    # Compute "chisq" - need to update to formal chisq
    mrqcof(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,y,sig,ndata,atry,ia,ma,covar,da,chisq,hextype,funcs)

    # if step succeeds value of chisq decreases: ochisq < chisq
    if chisq < ochisq:
        # decrease alamda by a factor of ten
        alamda *= 0.1
        # update chisq with the new, reduced, value
        ochisq = chisq 
        for j in range(mfit):
            for k in range(mfit):
                alpha[j][k] = covar[j][k]
            beta[j] = da[j]
        for l in range(ma):
            a[l] = atry[l]

    # else step does not succeed and chisq increases
    else:
        # Increase alamda by a factor of ten
        # and then loop back and try again
        alamda *= 10.0
        chisq = ochisq



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



# -------------------------------------------------------------------
# ---------------------------- MRQCOF -------------------------------
# -------------------------------------------------------------------
# ------------------------- Description from: -----------------------
# -------- Numerical Recipes: the art of scientific computing -------
# --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling --
# -------------------------------------------------------------------

# Subroutine MRQCOF(X, Y, SIG, NDATA, A, MA, LISTA, MFIT, ALPHA, BETA,
# NALP, CHISQ, FUNCS)
# Used by MRQMIN to evaluate the linearized fitting matrix ALPHA, and 
# vector BETA from (14.4.8) 

# called by mrqmin each iteration. Computes "chisq"
def mrqcof(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,y,sig,ndata,a,ia,ma,alpha,beta,chisq,hextype,funcs):

    mfit = 0
    dyda = numpy.zeros((ndata,ma))
    ymod = numpy.zeros(ndata)
    for j in range(ma):
        if ia[j] != 0:
            mfit += 1
    for j in range(mfit):
        for k in range(j):
            alpha[j][k] = 0.0
        beta[j] = 0.0

    chisq = 0.0
    funcs(d,r,itab,ltab,mtab,ntab,irk,d1,d2,d3,rho,shape,freqmin,ndata,a,ymod,dyda,ma,hextype)
    for i in range(ndata):
        sig2i = sig[i] * sig[i]
        dy = y[i] - ymod[i]
        j = -1
        for l in range(ma):
            if ia[l] != 0:
	        wt = dyda[l][i] * sig2i
                j += 1
                k = -1
                for m in range(l):
                    if ia[m] != 0:
                        k += 1
                        alpha[j][k] += wt * dyda[m][i]
	        beta[j] += dy * wt
 
        chisq += dy * dy * sig2i
  
    # chisq prints from here
    print('chisq={}'.format(100.0 * chisq))
    for j in range(1,mfit):
        for k in range(j-1):
            alpha[k][j] = alpha[j][k]


#/* ------------------------------------------------------------------- */
#/* ---------------------------- GAUSSJ ------------------------------- */
#/* ------------------------------------------------------------------- */
#/* ------------------------- Description from: ----------------------- */
#/* -------- Numerical Recipes: the art of scientific computing ------- */
#/* --W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T, Vetterling -- */
#/* ------------------------------------------------------------------- */

#/* Subroutine GAUSSJ(A, N, NP, B, M, MP)
# * Linear equation solution by Gauss-Jordan elimination, equation (2.1.1)
# * above. A is an input matrix of N by N elements, stored in an array of
# * physical dimensions NP by NP. B is an input matrix of N by N containing
# * the M right-hand side vectors, stored in an array of physical 
# * dimensions Np by MP. On output, A is replaced by its matrix inverse,
# * and B is replaced by the corresponding set of solution vectors. */
def gaussj(a, n, b, m):
    indxc = []
    indxr = []
    ipiv  = []
    for j in range(n):
        ipiv.append(0)
    for i in range(n):
        big = 0.0
        for j in range(n):
            if ipiv[j] != 1:
                for k in range(n):
                    if ipiv[k] == 0:
                        if math.fabs(a[j][k]) >= big:
                            big = math.fabs(a[j][k])
                            irow = j
                            icol = k
                    elif ipiv[k] > 1:
                        raise RuntimeError('GAUSSJ: Singular Matrix-1')
        ++(ipiv[icol]);
        if irow != icol:
            for l in range(n):
                a[irow][l], a[icol][l] = a[icol][l], a[irow][l]
            for l in range(m):
                b[irow][l], b[icol][l] = b[icol][l], b[irow][l]
        indxr.append(irow)
        indxc.append(icol)
        if a[icol][icol] == 0.0:
            raise RuntimeError('GAUSSJ: Singular Matrix-2')
        pivinv = 1.0 / a[icol][icol]
        a[icol][icol] = 1.0
        for l in range(n):
            a[icol][l] *= pivinv
        for l in range(m):
            b[icol][l] *= pivinv
        for ll in range(n):
            if ll != icol:
                dum = a[ll][icol];
                a[ll][icol] = 0.0
                for l in range(n):
                    a[ll][l] -= a[icol][l] * dum
                for l in range(m):
                    b[ll][l] -= b[icol][l] * dum
    for l in range(n-1,-1,-1):
        if indxr[l] != indxc[l]:
            for k in range(n):
                a[k][indxr[l]], a[k][indxc[l]] = a[k][indxc[l]], a[k][indxr[l]]

