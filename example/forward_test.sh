#! /bin/sh

set -x
# --------------------------------------------------------------------------

    d=8
    d1=4.44
    d2=4.44
    d3=6.44
   
    rho=2.72
    
    shape=1
    nfreq=20
    
    c44=113.0
    c11=26.0
 
    hextype=0
    
    rus_forward \
      d=$d d1=$d1 d2=$d2 d3=$d3 ns=2 nfreq=$nfreq\
        c11=$c11 c44=$c44 rho=$rho shape=$shape outeigen=0 eigenfile=eigenfunct\



