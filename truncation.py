import storage as S
import fit as F
import lmfit as lm
import numpy as np
import matplotlib.pyplot as plt
from sky_detect import chunks, clipped_stats
from scipy.optimize import leastsq
from scipy.ndimage.filters import median_filter as filt
from scipy import interpolate
import bisect

def translate_x(x, point, dupl_pos=-1, round_method='right'):
	"""converts a point in ordered list x to the appropriate index of x. 
	For duplicates, it returns the last index in the set of duplicates by default."""
	if point < x[0]: return 0
	if point > x[-1]: return len(x)
	if round_method is 'left': #finds closest by rounding down 
		i = bisect.bisect_right(x, point) - 1
	elif round_method is 'right': #finds closest by rounding up 
		i = bisect.bisect_right(x, point)
	return i

def overlap_chunks(l, n):
	return [l[i:i+n] for i in range(0, len(l)-n)]

def find_bulge_cutoff(profile, infoDF, fit_result, tol=0.2):
	"""returns end of bulge position and difference between D and B mags"""
	bulge, disc = F.sersic(fit_result.params, infoDF.zp, profile.R.values, comp=True, show=True)
	total = F.convert_mag(bulge + disc, infoDF.zp)
	disc = F.convert_mag(disc, infoDF.zp)
	bulge = F.convert_mag(bulge, infoDF.zp)
	diff = np.where((disc - total) < tol)[0]
	return diff[0], bulge[diff[0]] - disc[diff[0]]

def find_sky_cutoff(profile, infoDF, sigma=3.):
	sky_cutoff = sigma * infoDF.sky_unc
	print sky_cutoff

def clean_profile(profile, window, spline_smooth=0, median_smooth_perc=0.05):
	sl = slice(window[0], window[1])
	mask = ~profile.M[sl].isnull().values
	y = target.M[sl][mask].values
	x = target.R[sl][mask].values
	tck = interpolate.splrep(x,y,s=spline_smooth)
	xnew = np.linspace(x[0],x[-1], len(x))
	ynew = filt(interpolate.splev(xnew,tck,der=0), len(x) * median_smooth_perc)
	return ynew, xnew


def get_local_derivative(profile, window, chunk_size=2):
	M, R = clean_profile(profile, window)
	M_list = overlap_chunks(M, chunk_size)
	R_list = overlap_chunks(R, chunk_size)
	chunk_list = map(np.vstack, zip(M_list, R_list))
	straight = lambda p, chnk: chnk[0] - p[0] - (p[1]*chnk[1])
	fit_func = lambda c: leastsq(straight, [1.,1.], args=(c))[0][1]
	gradients = map(fit_func, chunk_list)
	med_Rs = map(np.median, R_list)
	return gradients, med_Rs

def find_break(profile, window, chunk_size=2):
	grad_list, R_list = get_local_derivative(profile, window, chunk_size)
	R_break = np.average(R_list, weights = 1. / (abs(np.array(grad_list))))
	return R_break, grad_list, R_list

def define_boundaries(profile, window, chunk_size=2, clip_sig=1.):
	R_break, gradients, med_Rs = find_break(profile, window, chunk_size)
	brk_pos = translate_x(med_Rs, R_break)
	inner, outer = gradients[:brk_pos], gradients[brk_pos:]
	means, stds = zip(*[clipped_stats(i, sig=clip_sig) for i in (inner, outer)])
	if means[0] >= means[1]: #is inner > outer?
		inner_bound = np.where(inner[::-1] > means[0] + stds[0])[0][0]
		outer_bound = np.where(outer < means[1] - stds[1])[0][0]
	else:
		inner_bound = np.where(inner[::-1] < means[0] - stds[0])[0][0]
		outer_bound = np.where(outer > means[1] + stds[1])[0][0]
	inner_bound = med_Rs[brk_pos-1::-1][inner_bound]
	outer_bound = med_Rs[brk_pos:][outer_bound]

	plt.plot(med_Rs, gradients)
	plt.axvline(R_break, linestyle='-')
	colours = ['r', 'g']
	for i in range(2):
		plt.axhline(y=means[i]+stds[i], linestyle='--', color=colours[i])
		plt.axhline(y=means[i]-stds[i], linestyle='--', color=colours[i])
		plt.axhline(y=means[i], linestyle='-', color=colours[i])
	plt.show()
	R_window = [profile.R.values[w] for w in window]
	return [R_window[0], inner_bound, outer_bound, R_window[1]], R_break, gradients, med_Rs

def fit_truncated(profile, infoDF, break_bounds, break_R):
	straight = lambda p, chnk: chnk[0] - p[0] - (p[1]*chnk[1])
	fit_func = lambda c: leastsq(straight, [1.,1.], args=(c))[0][1]
	gradients = map(fit_func, chunk_list)
	
if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 242
	target, info = tables[N], header.loc[N]
	result = F.fit_bulge_disc(target, info)
	lm.report_fit(result.params, show_correl=False)
	pos, delta = find_bulge_cutoff(target, info, result)

	bounds, brk, grads, grad_Rs = define_boundaries(target, [pos, -1])
	print bounds

	fig, ax, ax2 = F.plot_basic(result, target, info)
	for i in [ax, ax2]:
		i.axvline(x=brk)
		for b in bounds:
			i.axvline(x=b, linestyle='--')
	plt.show()
