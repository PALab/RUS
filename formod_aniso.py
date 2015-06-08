import argparse
import sys
import scipy
import scipy.linalg.lapack as lapack
from math import sqrt

def gamma_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, cm, shape, irk):
	gamma = [scipy.zeros((irk[i],irk[i])) for i in range(8)]
	c = stiffness(cm)
	for k in range(8):
		irs = 0
		for ik in range(k):
			irs += irk[ik]
		irf = irs + irk[k]

		irv = 0
		ir1 = irs
		while ir1 < irf:
			irh = 0
			ir2 = irs
			while ir2 < irf:
				i1 = itab[ir1]
				i2 = itab[ir2]
				l1 = ltab[ir1]
				l2 = ltab[ir2]
				m1 = mtab[ir1]
				m2 = mtab[ir2]
				n1 = ntab[ir1]
				n2 = ntab[ir2]
				gamma[k][irv][irh] = 0.0;
				if l1 > 0:
					j1 = 0
					if l2 > 0:
						j2 = 0
						l = l1 + l2 - 2
						m = m1 + m2
						n = n1 + n2
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * l1 * l2 * v
					if m2 > 0:
						j2 = 1
						l = l1 + l2 - 1
						m = m1 + m2 - 1
						n = n1 + n2
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * l1 * m2 * v
					if n2 > 0:
						j2 = 2
						l = l1 + l2 - 1
						m = m1 + m2
						n = n1 + n2 - 1
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * l1 * n2 * v
				if m1 > 0:
					j1 = 1
					if l2 > 0:
						j2 = 0
						l = l1 + l2 - 1
						m = m1 + m2 - 1
						n = n1 + n2
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * m1 * l2 * v
					if m2 > 0:
						j2 = 1
						l = l1 + l2
						m = m1 + m2 - 2
						n = n1 + n2
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * m1 * m2 * v
					if n2 > 0:
						j2 = 2
						l = l1 + l2
						m = m1 + m2 - 1
						n = n1 + n2 - 1
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * m1 * n2 * v
				if n1 > 0:
					j1 = 2
					if l2 > 0:
						j2 = 0
						l = l1 + l2 - 1
						m = m1 + m2
						n = n1 + n2 - 1
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * n1 * l2 * v
					if m2 > 0:
						j2 = 1
						l = l1 + l2
						m = m1 + m2 - 1
						n = n1 + n2 - 1
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * n1 * m2 * v
					if n2 > 0:
						j2 = 2
						l = l1 + l2
						m = m1 + m2
						n = n1 + n2 - 2
						v = volintegral(d1, d2, d3, l, m, n, shape)
						gamma[k][irv][irh] += c[i1][j1][i2][j2] * n1 * n2 * v
				ir2 += 1
				irh += 1
			ir1 += 1
			irv += 1
	return gamma

def stiffness(cm):
	c = scipy.zeros((3,3,3,3))
	for i in range(3):
		for j in range(3):
			if i == 0 and j == 0:
				a = 0
			elif i == 1 and j == 1:
				a = 1
			elif i == 2 and j == 2:
				a = 2
			elif (i == 1 and j == 2) or (i == 2 and j == 1):
				a = 3
			elif (i == 0 and j == 2) or (i == 2 and j == 0):
				a = 4
			else:
				a = 5
			for k in range(3):
				for l in range(3):
					if k == 0 and l == 0:
						b = 0
					elif k == 1 and l == 1:
						b = 1
					elif k == 2 and l == 2:
						b = 2
					elif (k == 1 and l == 2) or (k == 2 and l == 1):
						b = 3
					elif (k == 0 and l == 2) or (k == 2 and l == 0):
						b = 4
					else:
						b = 5
					c[i][j][k][l] = cm[a][b]
	return c

def e_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, rho, shape, irk):
	e = [scipy.zeros((irk[i],irk[i])) for i in range(8)]
	for k in range(8):
		irs = 0
		for ik in range(k):
			irs += irk[ik]
		irf = irs + irk[k]

		irv = 0
		ir1 = irs
		while ir1 < irf:
			irh = 0
			ir2 = irs
			while ir2 < irf:
				i1 = itab[ir1]
				i2 = itab[ir2]
				l1 = ltab[ir1]
				l2 = ltab[ir2]
				m1 = mtab[ir1]
				m2 = mtab[ir2]
				n1 = ntab[ir1]
				n2 = ntab[ir2]

				l = l1 + l2
				m = m1 + m2
				n = n1 + n2
				if i1 != i2:
					e[k][irv][irh] = 0.0
				else:
					e[k][irv][irh] = rho * volintegral(d1, d2, d3, l, m, n, shape);
				ir2 += 1
				irh += 1
			ir1 += 1
			irv += 1
	return e

def volintegral(d1, d2, d3, l, m, n, shape):
	if l % 2 == 1 or m % 2 == 1 or n % 2 == 1:
		return 0.0
	# ell. cylinder shape
	if shape == 1:
		return 4.0 * scipy.pi * d1**(l+1) * d2**(m+1) * d3**(n+1) / (n+1) * doublefact(l-1) * doublefact(m-1) / doublefact(l + m + 2)
	# spheroid shape
	if shape == 2:
		return 4.0 * scipy.pi * d1**(l+1) * d2**(m+1) * d3**(n+1) * doublefact(l-1) * doublefact(m-1) * doublefact(n-1) / doublefact(l + m + n + 3)
	# rp shape
	return 8.0 / ((l+1) * (m+1) * (n+1)) * d1**(l+1) * d2**(m+1) * d3**(n+1)

def doublefact(n):
	if n == -1:
		return 1
	elif n == 0:
		return 1
	elif n == 1:
		return 1
	else:
		return n * doublefact(n-2)

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
	help='type of hexagonal symmetry: 1=VTT, 2=HTI (only required if ns=5)')
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

NSTACK = 50
NSMALL = 7
FM = 7875
FA = 211
FC = 1663

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

itab = [0 for i in range(int(r))]
ltab = [0 for i in range(int(r))]
mtab = [0 for i in range(int(r))]
ntab = [0 for i in range(int(r))]

irk = [0 for i in range(8)]

ir = 0

# k == 0
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 0 and m % 2 == 0 and n % 2 == 0:
						itab[ir] = i
						ltab[ir] = l
						mtab[ir] = m
						ntab[ir] = n
						ir += 1
						irk[0] += 1
				elif i == 1:
					if l % 2 == 1 and m % 2 == 1 and n % 2 == 0:
						itab[ir] = i
						ltab[ir] = l
						mtab[ir] = m
						ntab[ir] = n
						ir += 1
						irk[0] += 1
				elif i == 2:
					if l % 2 == 1 and m % 2 == 0 and n % 2 == 1:
						itab[ir] = i
						ltab[ir] = l
						mtab[ir] = m
						ntab[ir] = n
						ir += 1
						irk[0] += 1
print("irk[0]=" + str(irk[0]))

# k == 1
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 0:
						if m % 2 == 0:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[1] += 1
				if i == 1:
					if l % 2 == 1:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[1] += 1
				if i == 2:
					if l % 2 == 1:
						if m % 2 == 0:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[1] += 1
print("irk[1]=" + str(irk[1]))

# k == 2
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[2] += 1
				if i == 1:
					if l % 2 == 1:
						if m % 2 == 0:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[2] += 1
				if i == 2:
					if l % 2 == 1:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[2] += 1
print("irk[2]=" + str(irk[2]))

# k == 3
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[3] += 1
				if i == 1:
					if l % 2 == 1:
						if m % 2 == 0:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[3] += 1
				if i == 2:
					if l % 2 == 1:
						if m % 2 == 1:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[3] += 1
print("irk[3]=" + str(irk[3]))

# k == 4
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 1:
						if m % 2 == 0:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[4] += 1
				if i == 1:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[4] += 1
				if i == 2:
					if l % 2 == 0:
						if m % 2 == 0:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[4] += 1
print("irk[4]=" + str(irk[4]))

# k == 5
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 1:
						if m % 2 == 0:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[5] += 1
				if i == 1:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[5] += 1
				if i == 2:
					if l % 2 == 0:
						if m % 2 == 0:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[5] += 1
print("irk[5]=" + str(irk[5]))

# k == 6
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 1:
						if m % 2 == 1:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[6] += 1
				if i == 1:
					if l % 2 == 0:
						if m % 2 == 0:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[6] += 1
				if i == 2:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[6] += 1
print("irk[6]=" + str(irk[6]))

# k == 7
for i in range(3):
	for l in range(d+1):
		for m in range(d-l+1):
			for n in range(d-l-m+1):
				if i == 0:
					if l % 2 == 1:
						if m % 2 == 1:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[7] += 1
				if i == 1:
					if l % 2 == 0:
						if m % 2 == 0:
							if n % 2 == 1:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[7] += 1
				if i == 2:
					if l % 2 == 0:
						if m % 2 == 1:
							if n % 2 == 0:
								itab[ir] = i
								ltab[ir] = l
								mtab[ir] = m
								ntab[ir] = n
								ir += 1
								irk[7] += 1
print("irk[7]=" + str(irk[7]))

e = e_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, args.rho, args.shape, irk)
gamma = gamma_fill(itab, ltab, mtab, ntab, r, d1, d2, d3, cm, args.shape, irk)
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
