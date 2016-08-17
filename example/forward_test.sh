#! /bin/sh

set -x
# --------------------------------------------------------------------------

    d=6
    d1=3.76
    d2=3.76
    d3=3.76
   
    rho=1.70
    
    shape=2
    nfreq=10
    
    c33=10.625
    c23=4
    c12=4.284
    c44=2.448
    c66=5.508
 
    hextype=1
    
    gprof /usr/local/bin/rus_forward \
      d=$d d1=$d1 d2=$d2 d3=$d3 hextype=$hextype ns=5 nfreq=$nfreq\
        c33=$c33 c23=$c23 c12=$c12 c44=$c44 c66=$c66 rho=$rho shape=$shape outeigen=0 eigenfile=eigenfunct\



