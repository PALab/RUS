import scipy

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
    if n ==  -1: return 1
    elif n == 0: return 1
    elif n == 1: return 1
    else: return n * doublefact(n-2)

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

def xindex(ax, x, index):
	# Copyright (c) Colorado School of Mines, 2011.
	# All rights reserved. 

	# Author:  Dave Hale, Colorado School of Mines, 12/25/89
	# Translated to Python: Paul Freeman, University of Auckland, 6/11/2015

	nx = len(ax)

	# initialize lower and upper indices and step
	lower = index
	if lower < 0:
		lower = 0
	if lower >= nx:
		lower = nx - 1
	upper = lower + 1
	step = 1

	# if x values increasing
	if ax[-1] > ax[0]:
		# find indices such that ax[lower] <= x < ax[upper]
		while lower > 0 and ax[lower] > x:
			upper = lower
			lower -= step
			step += step
		if lower < 0:
			lower = 0;
		while upper < nx and ax[upper] <= x:
			lower = upper
			upper += step
			step += step
		if upper > nx:
			upper = nx

		# find index via bisection
		middle = (lower + upper) // 2
		while middle != lower:
			if x >= ax[middle]:
				lower = middle
			else:
				upper = middle
			middle = (lower + upper) // 2

	# else, if not increasing
	else:
		# find indices such that ax[lower] >= x > ax[upper]
		while lower > 0 and ax[lower] < x:
			upper = lower
			lower -= step
			step += step
		if lower < 0:
			lower = 0;
		while upper < nx and ax[upper] >= x:
			lower = upper
			upper += step
			step += step
		if upper > nx:
			upper = nx

		# find index via bisection
		middle = (lower + upper) // 2
		while middle != lower:
			if x <= ax[middle]:
				lower = middle
			else:
				upper = middle
			middle = (lower + upper) // 2

	return lower



def index_relationship(d, r):
	itab = [0 for i in range(int(r))]
	ltab = [0 for i in range(int(r))]
	mtab = [0 for i in range(int(r))]
	ntab = [0 for i in range(int(r))]
	irk  = [0 for i in range(8)]

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

	return itab, ltab, mtab, ntab, irk
