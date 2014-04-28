import pandas as pd
import matplotlib.pyplot as plt
from astroML.plotting import hist
import numpy as np
from conf_test import wilson
# from matplotlib import rc
# rc('font',**{'family':'serif','serif':['Times New Roman']},)

def graph_binned(totalDF, propy, ax, region_list=None, sum_cond=None):
	t = totalDF[(totalDF.h2 > 0) & (totalDF.M1 < 50)]
	if region_list is None:
		region_list = [35., 68.]
	region_list = [0.] + region_list + [t.r_clust.max()]
	groups = []
	for i, v in enumerate(region_list[1:]):
		groups.append(t[(region_list[i] < t.r_clust) & (t.r_clust < v)])
	if sum_cond is not None:
		y, tots = zip(*[(len(v[v[propy] == sum_cond]), float(len(v))) for v in groups])
		yerr = zip(*[wilson(v, tots[i], 0.68) for i, v in enumerate(y)])
	else:
		y = [v[propy].mean() for v in groups]
	x = [r.r_clust.mean() for r in groups]
	ax.errorbar(x, np.array(y)/np.array(tots), yerr=yerr, fmt='-o')
	print region_list, yerr, np.array(y)/np.array(tots)

def get_bins(N, x1, x2):
	x = x1.append(x2)
	return np.linspace(x.min(), x.max(), N)

def histogram(totalDF, column, axis, binN=10, incl_untrunc=True):
	trunc = totalDF[(totalDF.untruncated == False) & (totalDF.h2 > 0) & (totalDF.M1 < 50)][column].dropna()
	if not incl_untrunc:
		untrunc = trunc
	else:
		untrunc = totalDF[(totalDF.untruncated == True) & (totalDF.h2 > 0) & (totalDF.M1 < 50)][column].dropna()
	bins = get_bins(binN, untrunc, trunc)
	wuntrunc, wtrunc = [1./len(untrunc)]*len(untrunc), [1./len(trunc)]*len(trunc)
	hist(trunc, bins, ax=axis, alpha=0.5, label='Truncated', weights=wtrunc)
	if incl_untrunc:
		hist(untrunc, bins, ax=axis, alpha=0.5, label='Untruncated', weights=wuntrunc)
	
def hists(totalDF, columns, size, binN, excl_untrunc_list):
	if type(binN) == int:
		binN = [binN] * len(columns)
	fig = plt.figure()
	replace = {'BD_ratio':'$B/D$', 'MB':'$\mu_e$[mag]', 'ReB':'$R_{e}$[arcsec]', 'nB':'$n$', 'M1':'$\mu_{0, inner}$[mag]', 'M2':'$\mu_{0, outer}$[mag]',\
	'h1':'$h_1$[arcsec]', 'h2':'$h_2$[arcsec]', 'brk_R':'$R_{brk}$[arcsec]', 'strength':'$S_h$', 'r_clust':'$r_{clust}$[arcsec]'}
	List = []
	for i, c in enumerate(columns):
		ax = fig.add_subplot(size[0], size[1], i+1)
		if c in excl_untrunc_list:
			histogram(totalDF, c, axis=ax, binN=binN[i], incl_untrunc=False)
		else:
			histogram(totalDF, c, axis=ax, binN=binN[i], incl_untrunc=True)
		ax.set_xlabel(replace[c], fontsize=17, labelpad=1.)
		ax.set_ylabel('fraction of S0s')
		List.append(ax)
	maxy = max(zip(*[i.get_ylim() for i in List])[1])
	for a in List:
		a.set_ylim([0., maxy])
	return fig, List

def type_V_clust(totalDF, ax):
	graph_binned(total, 'untruncated', ax, sum_cond=True)



if __name__ == '__main__':
	import storage as S
	data, header = S.import_directory()
	store_trunc = pd.HDFStore('fixed_truncations.h5')
	total = store_trunc['mainDF']

	# f, axs = hists(total, ['BD_ratio', 'ReB', 'nB', 'strength', 'r_clust', 'brk_R', 'M1', 'M2'], (2,4), 10, ['M1', 'M2', 'strength', 'brk_R'])
	# leg = axs[4].legend(loc='best', prop={'size':16}, fancybox=True)
	# leg.get_frame().set_alpha(0.5)
	# f.set_facecolor('white')
	# f.subplots_adjust(wspace=0.2, hspace=0.2)
	print total[(total.untruncated == True) & (total.h2 > 0) & (total.M1 < 50)][['BD_ratio', 'ReB', 'nB', 'r_clust', 'brk_R', 'M1', 'M2']].describe()

	# f = plt.figure()
	# f.set_facecolor('white')
	# ax = f.add_subplot(111)
	# graph_binned(total, 'untruncated', ax, sum_cond=False, region_list=[35., 68.])
	# ax.set_xlabel('$R_{cluster}$ [Mpc]')
	# ax.set_ylabel('local fraction of truncated S0s')
	
	# hist((header.sky_level - header.sky), ax=ax, bins='knuth', color='0.7')
	# ax.set_xlabel('$\Delta$ sky [$\gamma$ cts]')
	# ax.set_ylabel('Number')

	plt.show()