import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import storage as S
import bootstrapping as B
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

def graph_par(totalDF, x, y, ax=None):
	if ax is None:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	xerr, yerr = [None] * 2
	if x+'_low' in totalDF.columns:
		xerr = (totalDF[x+'_low'],totalDF[x+'_high'])
	if y+'_low' in totalDF.columns:
		yerr = (totalDF[y+'_low'], totalDF[y+'_high'])
	ax.errorbar(totalDF[x], totalDF[y], yerr=yerr, xerr=xerr, fmt='b.')
	ax.set_ylabel(y)
	ax.set_xlabel(x)
	return fig, ax


if __name__ == '__main__':
	# clust_data = pd.read_table(DIR+'\\fullcas_sky.dat', skiprows=0,skipinitialspace=True, escapechar='#',delimiter='\t')
	# Rs = clust_data.set_index('ID').r_clust
	# tables, header = S.import_directory()
	# store_trunc = pd.HDFStore('store_trunc_boot200.h5', 'r')
	store_fits = pd.HDFStore('store_again.h5', 'r')
	# boot_tr = store_trunc.trunc_boot
	# boot_ft = store_fits.bootstraps
	# boot_comb_tr = raw_combine(boot_tr, header)
	# boot_comb_ft = raw_combine(boot_ft, header)

	# strengths = get_strengths(boot_comb_tr)
	# truncation = classify_type(strengths)
	# fits = pd.DataFrame(get_stats_dict(boot_comb_ft)).join(Rs)

	# total = fits.join(truncation)
	total = store_fits['mainDF']
	
	graph_par(total, 'r_clust', 'brk_R')	
	plt.show()


	# print C.downbended.sum(), C.untruncated.sum(), C.upbended.sum()
	# C.upbended.hist()
	# plt.show()
	
	
