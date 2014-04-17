import numpy as np
import storage as S
import fit as F
import pandas as pd
import lmfit as lm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
rc('font',**{'family':'serif','serif':['Times New Roman']})

params={'axes.labelsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'figure.figsize':(8,6), 'figure.facecolor':'white'}
plt.rcParams.update(params)

def make_axes(N=4, norm_res=True, hist=True, MR_RATIO = 6, edge_pad = 0.05, level_pad = 0.02, sep = 0.):
	f = plt.figure()
	if N == 1:
		left = [0.+edge_pad] 
		right = [1.-edge_pad] 
		top = [1. - edge_pad]
		bottom = [0. + edge_pad]
		xvis =[True]
		yvis =[True]
	else:
		left = [0.+edge_pad, 0.50] *2
		right = [0.50, 1.-edge_pad] *2
		top = [1.-edge_pad]*2 + [0.50 - (level_pad/2)]*2
		if N > 2:
			bottom = [0.50+(level_pad/2)]*2 + [0.+edge_pad]*2
		else:
			bottom = [0.+edge_pad]*2
	hspace = [sep]*4
	yvis = [True, False]*4
	if N <=2: 
		xvis = [True]*2
	else: 
		xvis = [False]*2 + [True]*2
	div_dir = ['left', 'right']*2

	mains, resids = [], []
	for i in range(N):
		gs = gridspec.GridSpec(MR_RATIO+1, 1)
		gs.update(left=left[i], right=right[i], top=top[i], bottom=bottom[i], hspace=hspace[i])
		m = f.add_subplot(gs[0:MR_RATIO-1, :])
		r = f.add_subplot(gs[MR_RATIO-1:, :])
		m.xaxis.set_visible(False)
		if not yvis[i]:
			r.yaxis.tick_right()
			m.yaxis.tick_right()
		if not xvis[i]:
			r.set_xticklabels(r.get_xticklabels(), visible=False)
		r.yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=4))

		mains.append(m)
		resids.append(r)
	return f, mains, resids

def to_lmfitP(DF, index):
	series = DF.loc[index]
	P = lm.Parameters()
	for i,v in series.iteritems():
		P.add(i,v)
	return P

def plot_data(data, index, ax, fmt='b.', ecolor='0.75', capthick=1.):
	ax.errorbar(data[index].R.values, data[index].M.values, yerr=[data[index].M_err_down, 
		data[index].M_err_up], fmt=fmt, ecolor=ecolor, capthick=capthick)
	ax.invert_yaxis()

def trunc_plotter(truncP, X, ax, overlap_perc=0.05):
	straight = lambda p, x, sec: p['mu_'+sec].value +  (1.086 * x / p['h_'+sec].value)
	get_overlapx = lambda r, h: np.sqrt(r * r * h * h / ((h * h) + (1.086**2.)))

	overlap = (X[-1] - X[0]) * 0.5 * overlap_perc
	overlap_in_x = [get_overlapx(overlap, h) for h in [truncP['h_inner'].value, truncP['h_outer'].value]]
	break_R = (truncP['mu_inner'].value - truncP['mu_outer'].value) / ((1.086/truncP['h_outer'].value) - (1.086/truncP['h_inner'].value))
	pltR1 = np.linspace(X[0], break_R+(overlap_in_x[0]), 1000)
	pltR2 = np.linspace(break_R-(overlap_in_x[1]), X[-1], 1000)
	ax.plot(pltR1, straight(truncP, pltR1, 'inner'), 'b-', linewidth=2)
	ax.plot(pltR2, straight(truncP, pltR2, 'outer'), 'b-', linewidth=2)
	ax.axvline(x=truncP['brk_R'].value)

def plot_fits(data, infoDF, fitDF, truncDF, index, ax, bulge=True, pre_disc=True, inner_disc=True, outer_disc=True, total=True, linefmt=None):
	if linefmt is None:
		linefmt = ['g:', 'r--', 'k--']
	pltR = np.linspace(0,data[index].R.values[-1], 1000)
	p = to_lmfitP(fitDF, index)
	bulge, disc = F.sersic(p, infoDF.loc[index].zp, pltR, comp=True, show=True)
	ax.plot(pltR, F.convert_mag(bulge, infoDF.loc[index].zp), linefmt[0])
	ax.plot(pltR, F.convert_mag(disc, infoDF.loc[index].zp), linefmt[1])
	ax.plot(pltR, F.convert_mag(disc+bulge, infoDF.loc[index].zp), linefmt[2])

	p = to_lmfitP(truncDF, index)
	print p
	trunc_plotter(p, pltR, ax)

	ax.set_ylim(35,15)



if __name__ == '__main__':
	tables, info = S.import_directory()
	store = pd.HDFStore('store_again_100.h5')
	fits, truncs, infos = store.fits, store.truncations, store.info

	fig, ms, rs = make_axes(4)
	plot_data(tables, 0, ms[0])
	plot_fits(tables, infos, fits, truncs, 0, ms[0])
	plt.show()
