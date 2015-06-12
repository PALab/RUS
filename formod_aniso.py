import argparse
import sys
import scipy
import scipy.linalg.lapack as lapack
from math import sqrt
import rus

parser = argparse.ArgumentParser(description='Forward Algorithm')
parser.add_argument(
    '--d',
    type=int,
    required=True,
    help='order of polynomials used to fit the eigenmodes')
parser.add_argument(
    '--d1',
    type=float,
    required=True,
    help='first dimension in cm')
parser.add_argument(
    '--d2',
    type=float,
    required=True,
    help='second dimension in cm')
parser.add_argument(
    '--d3',
    type=float,
    required=True,
    help='third dimension in cm')
parser.add_argument(
    '--rho',
    type=float,
    required=True,
    help='density in g/cm3')
parser.add_argument(
    '--ns',
    type=int,
    required=True,
    choices=[2,3,5,6,9],
    help='number of stiffness coefficient given (determine the symmetry)')
parser.add_argument(
    '--hextype',
    nargs=1,
    type=int,
    choices=[1,2],
    help='type of hexagonal symmetry: 1=VTI, 2=HTI (only required if ns=5)')
parser.add_argument('--c11', nargs=1, type=float)
parser.add_argument('--c12', nargs=1, type=float)
parser.add_argument('--c13', nargs=1, type=float)
parser.add_argument('--c22', nargs=1, type=float)
parser.add_argument('--c23', nargs=1, type=float)
parser.add_argument('--c33', nargs=1, type=float)
parser.add_argument('--c44', nargs=1, type=float)
parser.add_argument('--c55', nargs=1, type=float)
parser.add_argument('--c66', nargs=1, type=float)
parser.add_argument(
    '--nfreq',
    nargs='?',
    const='10',
    default='10',
    type=int,
    help='gives you the nfreqth first eigen frequencies on the screen')
parser.add_argument(
    '--outeigen',
    nargs='?',
    const=sys.stdout,
    default=sys.stdout,
    type=argparse.FileType('w'),
    help='file where to put eigenvectors in double format')
parser.add_argument(
    '--shape',
    nargs='?',
    const='1',
    default='1',
    type=int,
    choices=[0,1,2],
    help='0=rectangle, 1=ellipsoidal cylinder, 2=spheroid')

args = parser.parse_args()

cm = [[0.0 for i in range(6)] for j in range(6)]

if args.ns == 2:
    # isotropic
    if args.c11 and args.c44:
        cm[0][0] = args.c11[0]
        cm[3][3] = args.c44[0]

        cm[0][0] = cm[0][0] / 100
        cm[0][1] = cm[0][1] / 100
        cm[3][3] = cm[3][3] / 100
        cm[1][1] = cm[2][2] = cm[0][0]
        cm[4][4] = cm[5][5] = cm[3][3]
        cm[0][2] = cm[1][2] = cm[0][1]
        cm[2][0] = cm[2][1] = cm[1][0] = cm[0][1]
    else:
        parser.error('please provide c11 and c44 values')
if args.ns == 3:
    # cubic
    if args.c11 and args.c12 and args.c44:
        cm[0][0] = args.c11[0]
        cm[0][1] = args.c12[0]
        cm[3][3] = args.c44[0]

        cm[0][0] = cm[0][0] / 100
        cm[0][1] = cm[0][1] / 100
        cm[3][3] = cm[3][3] / 100
        cm[1][1] = cm[2][2] = cm[0][0]
        cm[4][4] = cm[5][5] = cm[3][3]
        cm[0][2] = cm[1][2] = cm[0][1]
        cm[2][0] = cm[2][1] = cm[1][0] = cm[0][1]
    else:
        parser.error('please provide c11, c12, and c44 values')
if args.ns == 5:
    # hexagonal
    if args.hextype and args.hextype[0] == 1:
        # VTI
        if args.c33 and args.c23 and args.c12 and args.c44 and args.c66:
            cm[2][2] = args.c33[0]
            cm[1][2] = args.c23[0]
            cm[0][1] = args.c12[0]
            cm[3][3] = args.c44[0]
            cm[5][5] = args.c66[0]

            cm[2][2] = cm[2][2] / 100
            cm[1][2] = cm[1][2] / 100
            cm[0][1] = cm[0][1] / 100
            cm[3][3] = cm[3][3] / 100
            cm[5][5] = cm[5][5] / 100
            cm[0][0] = cm[1][1] = 2.0 * cm[5][5] + cm[0][1]
            cm[0][2] = cm[2][0] = cm[2][1] = cm[1][2]
            cm[1][0] = cm[0][1]
            cm[4][4] = cm[3][3]
        else:
            parser.error('please provide c33, c23, c12, c44, and c66 values')
    elif args.hextype and args.hextype[0] == 2:
        # HTI
        if args.c11 and args.c33 and args.c12 and args.c44 and args.c66:
            cm[0][0] = args.c11[0]
            cm[2][2] = args.c33[0]
            cm[0][1] = args.c12[0]
            cm[3][3] = args.c44[0]
            cm[5][5] = args.c66[0]

            cm[0][0] = cm[0][0] / 100
            cm[2][2] = cm[2][2] / 100
            cm[0][1] = cm[0][1] / 100
            cm[3][3] = cm[3][3] / 100
            cm[5][5] = cm[5][5] / 100
            cm[1][2] = cm[2][1] = cm[2][2] - 2.0 * cm[3][3]
            cm[0][2] = cm[1][0] = cm[2][0] = cm[0][1]
            cm[1][1] = cm[2][2]
            cm[4][4] = cm[5][5]
        else:
            parser.error('please provide c11, c33, c12, c44, and c66 values')
    else:
        parser.error('please provide hextype')
if args.ns == 6:
    # tetragonal
    if args.c11 and args.c33 and args.c23 and args.c12 and args.c44 and args.c66:
        cm[0][0] = args.c11[0]
        cm[2][2] = args.c33[0]
        cm[1][2] = args.c23[0]
        cm[0][1] = args.c12[0]
        cm[3][3] = args.c44[0]
        cm[5][5] = args.c66[0]

        cm[0][0] = cm[0][0] / 100
        cm[2][2] = cm[2][2] / 100
        cm[1][2] = cm[1][2] / 100
        cm[3][3] = cm[3][3] / 100
        cm[0][1] = cm[0][1] / 100
        cm[5][5] = cm[5][5] / 100
        cm[1][1] = cm[0][0]
        cm[0][2] = cm[2][0] = cm[1][2]
        cm[1][0] = cm[0][1]
        cm[2][1] = cm[1][2]
        cm[4][4] = cm[3][3]
    else:
        parser.error('please provide c11, c33, c23, c12, c44, and c66 values')
if args.ns == 9:
    # orthorhombic
    if args.c11 and args.c22 and args.c33 and args.c23 and args.c13 and args.c12 and args.c44 and args.c55 and args.c66:
        cm[0][0] = args.c11[0]
        cm[1][1] = args.c22[0]
        cm[2][2] = args.c33[0]
        cm[1][2] = args.c23[0]
        cm[0][2] = args.c13[0]
        cm[0][1] = args.c12[0]
        cm[3][3] = args.c44[0]
        cm[4][4] = args.c55[0]
        cm[5][5] = args.c66[0]

        cm[0][0] = cm[0][0] / 100
        cm[1][1] = cm[1][1] / 100
        cm[2][2] = cm[2][2] / 100
        cm[1][2] = cm[1][2] / 100
        cm[0][2] = cm[0][2] / 100
        cm[0][1] = cm[0][1] / 100
        cm[3][3] = cm[3][3] / 100
        cm[4][4] = cm[4][4] / 100
        cm[5][5] = cm[5][5] / 100
        cm[2][0] = cm[0][2]
        cm[1][0] = cm[0][1]
        cm[2][1] = cm[1][2]
    else:
        parser.error('please provide c11, c22, c33, c23, c13, c12, c44, c55, and c66 values')


# dimension of the problem
d = args.d
r = 3 * (d + 1) * (d + 2) * (d + 3) / 6
# half sample dimensions are used in calculations
d1 = args.d1 / 2.0
d2 = args.d2 / 2.0
d3 = args.d3 / 2.0

itab, ltab, mtab, ntab, irk = rus.index_relationship(d, r)
e = rus.e_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, args.rho, args.shape, irk)
gamma = rus.gamma_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, cm, args.shape, irk)
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
