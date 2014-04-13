import fit as F
import lmfit as lm
import numpy as np
import storage as S
import pandas as pd
import sys


def fit_basic(P, R, I, W, Z):
	fitter = lm.Minimizer(F.sersic, P, fcn_args=(Z, R, I, W))
	fitter.leastsq()
	return fitter

def bootstrap(P, profile, infoDF, size=1000, load_bar=False):
	resamples = np.random.randint(0, len(profile.R.values), (1000, len(profile.R.values)))
	cols = [name for name, par in P.iteritems()]
	DF = pd.DataFrame(None, None, cols, dtype=float)
	for i in range(size):
		mask = resamples[i]
		fit = fit_basic(P.copy(), profile.R.values[mask], profile.I.values[mask], profile.I_err.values[mask], infoDF.zp)
		d = {name: par.value for name, par in P.iteritems()}
		DF = DF.append(d, ignore_index=True)
		if load_bar:
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/size), (100*i/size)))
			sys.stdout.flush()
	return DF

if __name__ == '__main__':
	tables, headers = S.import_directory()
	N = 2
	target, info = tables[N], headers.loc[N]
	result = F.fit_bulge_disc(target, info)
	import time
	start = time.clock()
	data = bootstrap(result.params, target, info, size=100, load_bar=False)
	print "\n will be finished in: %.2fs" % ((time.clock() - start) * len(tables))

	hist, bins = np.histogram(data.MB.values, bins=50)
	width = 1. * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.bar(center, hist, align='center', width=width)
	plt.show()
