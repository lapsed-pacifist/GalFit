import fit as F
import lmfit as lm
import numpy as np
import storage as S
import pandas as pd
import sys
from scipy.optimize import leastsq
import time
import truncation as T

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

def sided_std(DF, mean):
	dev2 = (DF - mean)**2.
	var = (1./(len(DF)-1)) * dev2.sum(axis=0)
	return np.sqrt(var)

def find_stats(DF):
	mean = DF.mean(axis=0)
	upper = DF[DF >= mean]
	lower = DF[DF <= mean]
	return [sided_std(lower, mean), sided_std(upper, mean)], mean

def average_boot(boot_panel):
	for i, b in boot_panel.iteritems():
		stds, mean = find_stats(b) #series
		if i == 0:
			meanDF = pd.DataFrame(mean).T
			stdupDF = pd.DataFrame(stds[1]).T
			stddownDF = pd.DataFrame(stds[0]).T
		else:
			meanDF = meanDF.append(mean, ignore_index=True)
			stdupDF = stdupDF.append(stds[1], ignore_index=True)
			stddownDF = stddownDF.append(stds[0], ignore_index=True)
	return meanDF, stddownDF, stdupDF

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

def change_sky(profile, infoS, sigma=1.):
	sky = infoS.sky_unc * sigma
	profile.I = (profile.i_cts + sky) / infoS.scale / infoS.scale
	profile.M = infoS.zp - (2.5 * np.log10(profile.I))
	up = profile.I + profile.I_err
	down = profile.I - profile.I_err
	Mdown = infoS.zp - (2.5 * np.log10(abs(up)))
	Mup = infoS.zp - (2.5 * np.log10(abs(down)))
	profile['M_err_down'] = abs(profile.M - Mdown)
	profile['M_err_up'] = abs(Mup - profile.M)
	return profile

def bootstrap_sky(profile, infoS, truncS, boot_size=50, LoadBar=None):
	r = np.linspace(-1., 1., boot_size)
	bounds = truncS[['b0', 'b1', 'b2', 'b3']].values
	for i, sky in enumerate(r):
		if LoadBar is not None: LoadBar.time()
		adj_prof = change_sky(profile.copy(), infoS, sigma=sky * infoS.sky_unc)
		try:
			pars, brk = T.fit_truncated(adj_prof, infoS, bounds, truncS.brk_R, fix_brk=False)
			row = np.append(np.ravel(pars), np.array([brk, sky]))
		except TypeError:
			row = np.ones((6)) * np.nan	
		columns=['M1', 'h1', 'M2', 'h2', 'brk_R', 'sky_perc']
		row = {columns[i]: row[i] for i in range(6)}	
		if i == 0:
			DF = pd.DataFrame(data=row, index=[0])
		else:
			DF = DF.append(row, ignore_index=True)
		if LoadBar is not None: LoadBar.progress('%.i/%.i' % (i, boot_size))
	return DF

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
	