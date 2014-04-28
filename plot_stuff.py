"""
make axes returns a 2 lists of residual axes and main axes
takes a DF of all parameters and plots them
returns a table of all parameters
"""
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman']})
import pandas as pd
import storage as S
import fit as F
import lmfit as lm
import bootstrapping as B

def to_lmfitP(DF, index):
	series = DF.loc[index]
	P = lm.Parameters()
	for i,v in series.iteritems():
		P.add(i,v)
	return P

def make_axes(f=None):
	if f is None:
		f = plt.figure()
	gs = gridspec.GridSpec(6,6)
	axes, resids = [], []
	for i in range(2):
		axes.append(f.add_subplot(gs[0:5, (i*3):(i*3)+3]))
		resids.append(f.add_subplot(gs[5:6, (i*3):(i*3)+3], sharex=axes[i]))
		resids[i].yaxis.set_major_locator(MaxNLocator(prune='upper', nbins=4))
	plt.subplots_adjust(hspace=0., wspace=0.)
	for a in [resids[1], axes[1]]:
		a.yaxis.tick_right()
	for a in axes:
		plt.setp(a.get_xticklabels(), visible=False)
	# resids[0].invert_xaxis()
	return f, axes, resids

def plot_data(data, index, ax, color='#595959', marker='o', ms=4, ecolor='0.75', capthick=1., invert=False):
	R, M = data[index].R, data[index].M
	ax.errorbar(R, M, yerr=[data[index].M_err_down, data[index].M_err_up], color=color, ecolor=ecolor, capthick=capthick, fmt=marker, markersize=ms)
	ax.invert_yaxis()
	lim = [min(R), max(R)]
	if invert:
		ax.set_xlim(lim[::-1])
	else:
		ax.set_xlim(lim)
	ylim = [35, 15]
	ax.set_ylim([35, 15])
	plot_outliers(R, M, ylim[::-1], ax, 0.1, nan_orient='up')
	
def trunc_model(M1, M2, h1, h2, R, r_brk, overlap=5.):
	straight = lambda mu0, h, x: mu0 +  (1.086 * x / h)
	lowR, highR = R[R < (r_brk + overlap)], R[R >= (r_brk - overlap)]
	t1, t2 = straight(M1, h1, lowR), straight(M2, h2, highR)
	return t1, t2, lowR, highR, np.append(t1, t2)

def plot_fits(data, infoDF, fitDF, truncDF, index, ax):
	trc = truncDF.loc[index]
	Z = infoDF.zp[index]

	pltR = np.linspace(0,data[index].R.values[-1], 1000)
	p = to_lmfitP(fitDF, index)
	bulge, disc = F.sersic(p, Z, pltR, comp=True, show=True)
	ax.plot(pltR, F.convert_mag(bulge, Z), 'g--', linewidth=2)
	ax.plot(pltR, F.convert_mag(disc, Z), 'r:', linewidth=2)

	ts_lap = trunc_model(trc.inner_M, trc.outer_M, trc.inner_h, trc.outer_h, pltR, trc.brk_R)
	ts = trunc_model(trc.inner_M, trc.outer_M, trc.inner_h, trc.outer_h, pltR, trc.brk_R, overlap=0)
	ax.plot(ts_lap[2], ts_lap[0], color='#3399FF', linestyle='--' , linewidth=2)
	ax.plot(ts_lap[3], ts_lap[1], color='#000099', linestyle='--', linewidth=2)

	tot = F.convert_I(ts[4], Z)+bulge
	ax.plot(pltR, F.convert_mag(tot, Z), 'k-', linewidth=1)
	ax.set_ylim(35,15)

	

def plot_outliers(X, Y, limits, ax, len_scale=0.5, nan_orient=None):
	length = len_scale * (limits[1] - limits[0])
	for i, x in enumerate(X):
		c = 0
		if (np.isnan(Y[i])) and (nan_orient == 'up'):
			xytext=(x, limits[1] - length)
			xy=(x, limits[1])
			color = 'r'
			c = 1
		elif (np.isnan(Y[i])) and (nan_orient == 'down'):
			xytext=(x, limits[0] + length)
			xy=(x, limits[0])
			color = 'r'
			c = 1
		elif Y[i] > limits[1]:
			xytext=(x, limits[1] - length)
			xy=(x, limits[1])
			color = 'k'
			c = 1
		elif Y[i] < limits[0]:
			xytext=(x, limits[0] + length)
			xy=(x, limits[0])
			color = 'k'
			c = 1
		if c == 1:
			ax.annotate("", xy=xy, xycoords='data', xytext=xytext, textcoords='data', arrowprops=dict(arrowstyle="simple", color=color))

def plot_residuals(data, infoDF, fitDF, truncDF, index, ax, limits=(-2, 2), norm=True, color='#595959', marker='o', ms=4, ecolor='0.75', capthick=1.):
	trc = truncDF.loc[index]
	Z = infoDF.zp[index]
	model = trunc_model(trc.inner_M, trc.outer_M, trc.inner_h, trc.outer_h, data[index].R.values, trc.brk_R, overlap=0)[-1]
	res = data[index].M - model
	err = data[index].M_err_up
	mask = data[index].M > 0
	redchi = np.sum((res[mask] / err[mask]) ** 2.) / (7 + len(res))
	s = '$\chi^{2}_{\\nu} = %.2f$' % redchi
	ax.text(0.98, 0.98, s, va='top', ha='right', transform=ax.transAxes)
	if norm:
		y = res / err
		ax.plot(data[index].R, y, 'b.')
	
	else:
		y = res
		ax.errorbar(data[index].R, y, yerr=err, color=color, ecolor=ecolor, capthick=capthick, fmt=marker, markersize=ms)
		
	ax.set_ylim(limits)

def draw_bounds(data, truncDF, infoDF, index, ax, res_ax):
	for i in range(3):
		b = truncDF['b'+str(i)][index]
		ax.axvline(x=b, linestyle=':', color='b', ymax=0.1)
		ax.axvline(x=b, linestyle=':', color='b', ymin=0.9)
		res_ax.axvline(x=b, linestyle=':', color='b', ymin=0.9)
	sky = data[index].R.values[-infoDF.sky_pos[index]]
	ax.axvline(x=sky, linestyle=':', color='b', ymax=0.1)
	ax.axvline(x=sky, linestyle=':', color='b', ymin=0.9)
	res_ax.axvline(x=sky, linestyle=':', color='b', ymin=0.9)

def draw_sky(infoDF, index, ax):
	sky = infoDF.mu_crit[index]
	ax.axhspan(sky, ax.get_ylim()[0], facecolor='r', alpha=0.1)

def join_text(mean, stds, cols='F'):
	if cols is 'F': 
		cols = ['BD_ratio', 'MB', 'ReB', 'nB', 'h']
		replace = [r'B/D', r'\mu_{e}', r'R_{e}', r'n', r'h']
	else:
		cols = ['M1', 'M2', 'h1', 'h2']
		replace = [r'\mu_{0, inner}', r'\mu_{0, outer}', r'h_{inner}', r'h_{outer}']

	str2_map = lambda x: '%.2f' % x
	str1_map = lambda x: '%.1f' % x
	m = mean.map(str2_map)
	u, d = [(i).map(str2_map) for i in stds]
	txt = m + r'^{+' + u + r'}_{-' + d + r'}$'
	txt.rename({a:b for a,b in zip(cols, replace)}, inplace=True)
	return [r'$'+n+r'='+v for n, v in txt.iteritems() if n in replace]



def add_text(fitDF, fit_bootP, truncDF, trunc_bootP, index, ax, invert=False):	
	fit_stds, fit_mean = B.find_stats(fit_bootP[index])
	fit_stds[0]['h'] = fit_stds[0].ReD / 1.678
	fit_stds[1]['h'] = fit_stds[1].ReD / 1.678
	fit_mean['h'] = fit_mean.ReD / 1.678
	trunc_stds, trunc_mean = B.find_stats(trunc_bootP[index]) # get stds
	ftxt = join_text(fit_mean, fit_stds)
	ttxt = join_text(trunc_mean, trunc_stds, 'T')
	s = '\n'.join((ftxt+ttxt))
	pr = dict(boxstyle='round', facecolor='white', alpha=0.7)
	if invert: 
		ax.text(0.02, 0.98, s, ha='left', va='top', transform=ax.transAxes, fontsize=14, bbox=pr)
	else:
		ax.text(0.98, 0.98, s, ha='right', va='top', transform=ax.transAxes, fontsize=14, bbox=pr)

	

def plot_one(table_list, infoDF, fitDF, truncDF, index, ax, rs_ax, invert=True):
	plot_data(table_list, index, ax, invert=invert)
	plot_fits(table_list, infoDF, fitDF, truncDF, index=index, ax=ax)
	plot_residuals(table_list, infoDF, fitDF, truncDF, index=index, ax=rs_ax)
	draw_bounds(table_list, truncDF, infoDF, index, ax, rs_ax)
	draw_sky(infoDF, index, ax)

def draw_camera(table_list, infoDF, fitDF, truncDF, fit_bootP, trunc_bootP, cam_name, ID, axs, rs_axs, fig):
	index = infoDF[(infoDF.ID == ID) & (infoDF.cam == cam_name)].index
	plot_one(table_list, infoDF, fitDF, truncDF, index[0], axs[0], rs_axs[0], invert=True)
	plot_one(table_list, infoDF, fitDF, truncDF, index[1], axs[1], rs_axs[1], invert=False)
	axs[0].set_ylabel('$\mu$ [mag arcsec$^{-1}$]')
	rs_axs[0].set_ylabel('$|\mu|$')
	rs_axs[0].set_xlabel('R [arcsec]', horizontalalignment='center', verticalalignment='top')
	ticklab = rs_axs[0].xaxis.get_ticklabels()[0]
	trans = ticklab.get_transform()
	rs_axs[0].xaxis.set_label_coords(1., -0.2, transform=trans)
	maxR = max([i.get_xlim()[1] for i in rs_axs])
	rs_axs[0].set_xlim([maxR, 0])
	rs_axs[1].set_xlim([0, maxR])
	fig.suptitle(str(ID)+' '+str(cam_name))
	add_text(fitDF, fit_bootP, truncDF, trunc_bootP, index[0], axs[0], invert=True)
	add_text(fitDF, fit_bootP, truncDF, trunc_bootP, index[1], axs[1], invert=False)


def convert_to_lmfit(row, columns):
	P = lm.Parameters()
	for c in columns:
		P.add(c, row[c])
	P.add('deltaRe', vary=True)
	return P


if __name__ == '__main__':
	tables, header = S.import_directory()
	store = pd.HDFStore('store_again.h5')
	store_trunc = pd.HDFStore('fixed_truncations.h5')
	trunc = store_trunc.truncations
	fit = store.fits
	main = store_trunc.mainDF.reset_index()
	m = main[main.upbended == True].iloc[36]
	print m

	cam_name = 'mega'
	ID = m['index']
	f = plt.figure()
	f.set_facecolor('white')
	f, axes, rs = make_axes(f)
	draw_camera(tables, header, fit, trunc, store.bootstraps, store_trunc.trunc_boot, str(cam_name), ID, axes, rs, f)
	plt.show()

	# for ID in header.ID.unique():
	# 	for cam_name in header.cam.unique():
	# 		f.clf()
	# 		f, axes, rs = make_axes(f)
	# 		draw_camera(tables, header, fit, trunc, str(cam_name), ID, axes, rs, f)
	# 		plt.draw()
