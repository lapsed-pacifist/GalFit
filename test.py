import numpy as np


def bin(x, bins=None, n=None):
	"""rebins data and produces mask for other arrays
	Bins is a list of values corresponding to boundaries of the bins (ends have to be included!)"""
	if bins is None:
		bins = np.arange(x[0], x[-1], n)

	chks = [np.where((x <= x[bins[i+1]]) | x >= x[bins[i]]) for i in range(len(bins)-1)]
	return chks



x = np.linspace(0, 20)
n = 5

print x
print bin(x, n=n)