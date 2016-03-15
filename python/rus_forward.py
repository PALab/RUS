import rus_parser
import sys
import scipy
import scipy.linalg.lapack as lapack
from math import sqrt
import rus_tools as rus

def start(args):

    order   = args.order
    density = args.rho
    shape   = args.shape

    # half sample dimensions are used in calculations
    dimension1 = args.d1 / 2.0
    dimension2 = args.d2 / 2.0
    dimension3 = args.d3 / 2.0
    dimensions = [dimension1, dimension2, dimension3]

    # general size of the problem
    a = order + 1
    b = order + 2
    c = order + 3
    problem_size = 3 * a * b * c / 6

    cm = rus.calc_forward_cm(args)

    tabs, irk = rus.index_relationship(order,problem_size)

    e = rus.e_fill(tabs,dimensions,density,shape,irk)

    gamma = rus.gamma_fill(tabs,dimensions,cm,shape,irk)

    print("done preparing matrices")

    if args.outeigen == None:
        jobz = 'N'
    else:
        jobz = 'V'

    w = []
    for k in range(8):
        # lapack routine
        a, w_temp, info = lapack.dsygv(gamma[k], e[k], itype=1, jobz=jobz, uplo='U');  
        w.append(w_temp)

    wsort = scipy.zeros(problem_size)
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
