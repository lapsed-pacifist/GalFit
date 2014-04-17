import fit as F
import lmfit as lm
import numpy as np
import storage as S
import pandas as pd
import sys
from scipy.optimize import leastsq
import time

MAXFEV = 1000

def copy_params(parameters, trans_vary=True):
	new = lm.Parameters()
	for p in parameters.values():
		new.add(p.name, value=p.value, min=p.min, max=p.max)
		if trans_vary:
			new[p.name].vary=p.vary
	return new

def fit_basic(P, R, I, W, Z):
	fitter = lm.Minimizer(F.sersic, P, fcn_args=(Z, R, I, W))
	fitter.leastsq(maxfev=MAXFEV)
	return fitter

def bootstrap(P, profile, infoDF, size=1000, load_bar=False):
	resamples = np.random.randint(0, len(profile.R.values), (size, len(profile.R.values)))
	cols = [name for name, par in P.iteritems()]
	DF = pd.DataFrame(None, None, cols, dtype=float)
	ev = []
	for i in range(size):
		mask = resamples[i]
		fit = fit_basic(P.copy(), profile.R.values[mask], profile.I.values[mask], profile.I_err.values[mask], infoDF.zp)
		evals = fit.nfev
		ev.append(evals)
		d = {name: par.value for name, par in P.iteritems()}
		d.update({'nfev':evals})
		DF = DF.append(d, ignore_index=True)
		if load_bar:
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/size), (100*i/size)))
			sys.stdout.write(" max nfev:%i" % max(ev))
			sys.stdout.flush()
	return DF

def bootstrap_line(P, R, M, MW, size=1000, load_bar=False, label='in'):
	"""P is a numpy array of 2 parameters mu, h for a straight line"""
	mask = ~M.isnull()
	print MW[mask]
	resamples = np.random.randint(0, len(R), (size, len(R)))
	cols = ['mu0', 'h']
	DF = pd.DataFrame(None, None, cols, dtype=float)
	straight = lambda p, x,y,w: (y - (p[0] +  (1.086 * x / p[1]))) / w
	for i in range(size):
		mask = resamples[i]
		fit = leastsq(straight, P.copy(), args=(R[mask].values,M[mask].values,MW[mask].values))[0]
		d = {cols[i]+'_'+str(label):fit[i] for i in range(2)}
		DF = DF.append(d, ignore_index=True)
		if load_bar:
			sys.stdout.write('\r')
			sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*i/size), (100*i/size)))
			sys.stdout.flush()
	return DF

def bootstrap_trunc(P_list, R, M, MW, brk_pos, size=1000, load_bar=False):
	"""P_list is a list of 2 arrays containing the fit parameters [mu0, h]"""
	inner = bootstrap_line(P_list[0], R[:brk_pos+1], M[:brk_pos+1], MW[:brk_pos+1], size, load_bar, 'in')
	outer = bootstrap_line(P_list[1], R[brk_pos:], M[brk_pos:], MW[brk_pos:], size, load_bar, 'out')
	return pd.concat((inner, outer), axis=1)

def sided_std(series, mean):
	dev2 = (series - mean)**2.
	var = (1./(len(series)-1)) * dev2.sum()
	return np.sqrt(var)

def find_stats(series):
	mean = series.mean()
	upper = series[series >= mean].values
	lower = series[series <= mean].values
	return [sided_std(lower, mean), mean, sided_std(upper, mean)]

# def dist(x, mean, std):
# 	pre = 1./ std / np.sqrt(2*np.pi)
# 	exponent = (x - mean) **2. / (2.* std * std)
# 	return pre * np.exp(-exponent)

# def sided_dist(mean, std_low, std_high, start=None, stop=None, num=1000):
# 	if start is None: start = std_low
# 	if stop is None: stop = std_high
# 	X_low = np.linspace(start, mean, num/2)
# 	X_high = np.linspace(mean, stop, num/2)
# 	low = dist(X_low, mean, std_low)
# 	high = dist(X_high, mean, std_high)
# 	return np.append(low, high[1:]), np.append(X_low, X_high[1:])

def group_lists(stat_list):
	"""takes any number of lists of [-std, mean +std] and groups them if they overlap"""
	stat_list.sort(key=lambda k: (k[1]))
	comb_list = []
	block_no = 0
	for i,s in enumerate(stat_list):
		if i == 0:
			comb_list.append([s]) 
		else:
			Block =	comb_list[block_no]
			bound = max(map(lambda x: x[1]+x[2], Block))
			if s[1] - s[0] <= bound:
				comb_list[block_no].append(s)
			else:
				comb_list.append([s])
				block_no += 1
	return comb_list


def combine_boot(*boot):
	stats = [find_stats(b).append(i) for i,b in enumerate(boot)]
	combinations = group_lists(stats)
	boot_list = []
	for c in combinations:
		boots = [boot[i[3]] for i in c]
		boot_list.append(np.hstack(tuple(boots)))
	return boot_list


if __name__ == '__main__':
	tables, headers = S.import_directory()
	N = 36
	target, info = tables[N], headers.loc[N]
	result = F.fit_bulge_disc(target, info)
	data = bootstrap(result.params, target, info, size=1000, load_bar=True)
	print data.describe()


	hist, bins = np.histogram(data.MB.values, bins=50)
	width = 1. * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.bar(center, hist, align='center', width=width)
	plt.show()
	