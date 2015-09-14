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

# global variable
_dqkpart_seed = 0

# send to parser
args = p.inverse_parser(sys.argv)

print('d={}'.format(args.order))
print('shape={}'.format(args.shape))
print('ns={}'.format(args.ns))
print('hextype={}'.format(args.hextype))
print('d1={0:.6f}'.format(args.d1))
print('d2={0:.6f}'.format(args.d2))
print('d3={0:.6f}'.format(args.d3))
print('rho={0:.6f}'.format(args.rho))
print('freqmin={0:.6f}'.format(args.freqmin))
print('freqmax={0:.6f}'.format(args.freqmax))

# print out initial guess
for k,v in args.a.iteritems():
    print('{0:.6f}'.format(v))
    args.a[k] = v / 100

d  = args.order
d1 = args.d1 / 2.0  # half sample dimensions are used in calculations
d2 = args.d2 / 2.0
d3 = args.d3 / 2.0 

# dimension of the problem
r = 3 * (d + 1) * (d + 2) * (d + 3) / 6

if args.input_file == 'sample/default_frequencies':
    print('!! no frequency file specified - sample file will be used !!')

# get measured frequencies from file
try:
    f = open(args.input_file, "rU")
except IOError:
    print('Could not open frequency file.')
    sys.exit(-1)
else:
    nfreq = int(f.readline())
    print('nfreq={}'.format(nfreq))
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
        print(' freq={0:.6f}'.format(freq[-1]))
    f.close()

if len(freq) != nfreq:
    print('Unexpected number of frequencies.')
    print('Expected {} but read {}'.format(nfreq, len(freq)))
    sys.exit(-1)

# relationship between ir and l,m,n - filling tables
tabs, irk = rus.index_relationship(d, r)

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
alamda = -1.0
orig_cxx_values = args.a.copy()
for i in range(args.iterations):
    rus.mrqmin(d,r,tabs,irk,d1,d2,d3,args.rho,args.shape,args.freqmin,y,sig,ndata,args.a,ia,args.ns,covar,alpha,args.hextype,alamda)
    print('iter #{}'.format(i))
    for k in sorted(args.a.keys()):
        print('{}'.format(100 * args.a[k])) # print estimated cxx values

print('\nThis calculation can be executed again with the following command:')
print('python {} --order {} --shape {} --ns {} --hextype {} --d1 {} --d2 {} --d3 {} --rho {} --freqmin {} --freqmax {} --iterations {} {}'.format(sys.argv[0], args.order, args.shape, args.ns, args.hextype, args.d1, args.d2, args.d3, args.rho, args.freqmin, args.freqmax, args.iterations, ' '.join(('--' + k + ' ' + str(v*100)) for k,v in orig_cxx_values.iteritems())))



