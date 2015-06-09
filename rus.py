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
