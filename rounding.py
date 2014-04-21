import numpy as np
from decimal import Decimal 

def round_to(x, sig=2):
	dec = sig-int(np.floor(np.log10(x)))-1
	a = np.round(x, dec)
	return str(a), dec

def round_err(value, err):
	re = round_to(err, 1)
	rv = str(np.round(float(v), re[1]))
	diff = len(re[0]) - len(rv)
	if diff > 0:
		rv += '0' * diff
	return rv, re[0]

def round_2err(value, err1, err2):
	first = round_err(value, float(err1))
	sec = round_err(value, float(err2))
	if float(sec[0]) > float(first[0]):
		return sec[0], first[1], sec[1]
	else:
		return first[0], first[1], sec[1]

if __name__ == '__main__':
	v, e1, e2 = 2, 3, 10
	print round_2err(v, e1, e2)