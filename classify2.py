import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import storage as S
import bootstrapping as B
from fit import get_b_n
from scipy.special import gammainc, gamma
from fit import convert_I
import conf_test
from  scipy.optimize import leastsq
from scipy.stats import pearsonr
DIR = 'repository'

def get_stats_dict(bootP):
	"""returns a DF with a std+- and value for each column in panel"""
	cols = bootP[bootP.items[0]].columns
	d = {}
	for c in cols:
		df = bootP.minor_xs(c)
		stds, mean = B.find_stats(df)
		d.update({c:mean, c+'_low':stds[0], c+'_high':stds[1]})
	return d

def get_strengths(bootP):
	"""returns a DF of mean strength with its std+-"""
	d = get_stats_dict(bootP)
	df_str = (bootP.minor_xs('h2') / bootP.minor_xs('h1')) -1.
	stds, mean = B.find_stats(df_str)
	d.update({'strength':mean, 'strength_low':stds[0], 'strength_high':stds[1]})
	return pd.DataFrame(d)

def index_by_ID(df, infoDF, cols=None):
	if cols is None:
		cols = ['ID','cam','ax']
	return pd.concat([df, header[cols]], axis=1).set_index(cols).sortlevel()
	
def raw_combine(bootP, infoDF):
	"""returns a DF of mean strength with its std+- regardless of compatibility"""
	inds = [infoDF[infoDF.ID == i].index for i in infoDF.ID.unique()]
	for i, group in enumerate(inds): # for each galaxy
		DF = bootP[group[0]]
		for j in group[1:]:
			DF = DF.append(bootP[j], ignore_index=True)
		if i == 0:
			wp = pd.Panel({infoDF.ID[group[0]]:DF})
		else:
			wp[infoDF.ID[group[0]]] = DF
	return wp

def classify_type(strengthDF):
	h1_std_ratio = (strengthDF.h1 + strengthDF.h1_high) / (strengthDF.h1 - strengthDF.h1_low)
	Sup, Sdown = strengthDF.strength + strengthDF.strength_high, strengthDF.strength - strengthDF.strength_low
	is_untrunc = ~(((Sup + 1) < (1 / h1_std_ratio)) | ((Sdown + 1) > h1_std_ratio))
	is_upbend = (strengthDF.strength >= 0) & (~is_untrunc)
	is_downbend = (strengthDF.strength <= 0) & (~is_untrunc)
	types = pd.DataFrame({'untruncated': is_untrunc, 'upbended':is_upbend, 'downbended':is_downbend})
	return strengthDF.join(types)

def bin(x, bins=None, n=None):
	"""rebins data and produces mask for other arrays
	Bins is a list of values corresponding to boundaries of the bins (ends have to be included!)"""
	if bins is None:
		bins = range(x[0], x[-1], n)

	chks = [np.where(x[bins[i]] <= x <= x[bins[i+1]]) for i in range(len(bins)-1)]

def correl(X, Y, Xerr, Yerr):
	straight = lambda p, x, y, w: (y - (p[0]* x) - p[1]) / w
	pars = leastsq(straight, [1.,1.], args=(X, Y, Yerr))[0]
	return pars[0]

def graph_par(totalDF, x, y, bins=None, ax=None, err=True, truncated_only=False, fmt='b.', P=True):
	if truncated_only:
		t = totalDF[totalDF.untruncated != False]
	else:
		t = totalDF
	if ax is None:
		f = plt.figure()
		ax = fig.add_subplot(111)
	else:
		f = None
	xerr, yerr = [None] * 2
	if x+'_low' in t.columns:
		xerr = (t[x+'_low'],t[x+'_high'])
	if y+'_low' in t.columns:
		yerr = (t[y+'_low'], t[y+'_high'])
	if err:
		ax.errorbar(t[x], t[y], yerr=yerr, xerr=xerr, fmt=fmt, ecolor='0.75')
	else:
		ax.plot(t[x], t[y], 'b.')
	replace = {'BD_ratio':'$B/D$', 'MB':'$\mu_e$[mag]', 'ReB':'$R_{e}$[arcsec]', 'nB':'$n$', 'M1':'$\mu_{0, inner}$[mag]', 'M2':'$\mu_{0, outer}$[mag]',\
	'h1':'$h_1$[arcsec]', 'h2':'$h_2$[arcsec]', 'brk_R':'$R_{brk}$[arcsec]', 'strength':'$S_h$', 'r_clust':'$r_{clust}$[arcsec]'}
	ax.set_ylabel(replace[y])
	ax.set_xlabel(replace[x])
	if fig is None:
		return ax
	else:
		return f, ax

def add_corr_coeff(DF, x, y, ax, text='', pos_y=0):
	mask = ~DF[x].isnull() & ~DF[y].isnull()
	X = DF[x][mask].values
	Y = DF[y][mask].values
	P = pearsonr(X, Y)
	ax.text(0.98, 0.98-pos_y, '$P_{%s} = %.2f$' % (text, P[0]), transform=ax.transAxes, ha='right', va='top')

def new_BD_ratio(P, zp):
	MB, ReB, nB, M1, M2, h1, h2, R_brk = P.values
	I01, I02, IeB = convert_I(M1, zp), convert_I(M2, zp), convert_I(MB, zp)
	disc1 = I01 * h1 * h1 * gammainc(2., R_brk / h1)
	disc2 = I02 * h2 * h2 * (1 - gammainc(2., R_brk / h2))
	bn = get_b_n(nB)
	bulge = np.exp(bn) * IeB * ReB * ReB * gamma((2 * nB) + 1) / (bn**(2.*nB))
	return bulge / (2* (disc2 + disc1))

def get_BD_ratios(totalS, zp):
	cols = ['MB', 'ReB', 'nB', 'M1', 'M2', 'h1', 'h2', 'R_brk']
	BD = new_BD_ratio(totalS[cols], zp)
	low = new_BD_ratio(totalS[[i+'_low' for i in cols]], zp)
	high = new_BD_ratio(totalS[[i+'_high' for i in cols]], zp)
	return BD, low, high

def graph_type_hist(totalDF):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	mine = totalDF[['untruncated', 'downbended', 'upbended']].sum().values
	erwin = [0.458*55., 0., (1-0.458)*55.] # unt, down, up
	erwin_errs = [conf_test.wilson(i, 55, 0.68) for i in erwin]
	mine_errs = [conf_test.wilson(i, float(len(totalDF)), 0.68) for i in mine]
	
	ind = np.arange(3)
	width = 0.5
	rects1 = ax.bar(ind, np.array(mine)/np.sum(mine), width, color='0.8', yerr=zip(*mine_errs))
	rects2 = ax.bar(ind+width, np.array(erwin)/55., width, color='0.3', yerr=zip(*erwin_errs))


if __name__ == '__main__':
	# clust_data = pd.read_table(DIR+'\\fullcas_sky.dat', skiprows=0,skipinitialspace=True, escapechar='#',delimiter='\t')
	# Rs = clust_data.set_index('ID').r_clust
	tables, header = S.import_directory()
	store_trunc = pd.HDFStore('fixed_truncations.h5')

	# store_fits = pd.HDFStore('store_again.h5', 'r')

	# boot_tr = store_trunc.trunc_boot
	# boot_ft = store_fits.bootstraps
	# boot_comb_tr = raw_combine(boot_tr, header)
	# boot_comb_ft = raw_combine(boot_ft, header)

	# strengths = get_strengths(boot_comb_tr)
	# truncation = classify_type(strengths)
	# fits = pd.DataFrame(get_stats_dict(boot_comb_ft)).join(Rs)

	# total = fits.join(truncation)
	total = store_trunc['mainDF']
	# print total.loc[1237667322723827758]
	# print total.BD_ratio.iloc[0]
	# print get_BD_ratios(total.iloc[0], header.zp[0])
	
	fig = plt.figure()	
	fig.set_facecolor('white')
	x, y = 'r_clust', 'h2'
	ax1 = fig.add_subplot(111)

	t = total[(total.h1 > 0) & (total.untruncated==True)]
	add_corr_coeff(t, x, y, ax1, text='type-I', pos_y=0)
	graph_par(t, x, y, ax=ax1, P=False)

	t = total[(total.h1 > 0) & (total.untruncated==False)]
	graph_par(t, x, y, ax=ax1, fmt='r.', P=True)
	add_corr_coeff(t, x, y, ax1, text='type-III', pos_y=0.05)

	t = total[(total.h1 > 0)]
	add_corr_coeff(t, x, y, ax1, text='total', pos_y=0.1)

	# t = total[(total.h1 > 0) & (total.untruncated==False)]
	# t[['MB','ReB', 'nB', 'ReD']].hist()
	# plt.title('truncated')

	# t = total[(total.h1 > 0) & (total.untruncated==True)]
	# t.MB.hist(ax=fig.add_subplot(122))
	# plt.title('untruncated')

	# ff = lambda x: '%.2f' % (x)
	# print t[['MB','ReB','nB', 'BD_ratio']].describe().to_latex(float_format=ff)
	# plt.ylabel('Number')
	# plt.xlabel('$h_{outer}$[arcsec]', labelpad=0.5)
	# graph_type_hist(total)
	plt.show()



