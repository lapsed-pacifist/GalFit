import storage as S
import fit as F
import lmfit as lm
import numpy as np
import matplotlib.pyplot as plt
from sky_detect import chunks
from scipy.optimize import leastsq
from scipy.ndimage.filters import median_filter as filt
from scipy import interpolate

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

def clean_profile(profile, window, spline_smooth=0, median_smooth_perc=0.1):
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
	gradients, med_Rs = get_local_derivative(profile, window, chunk_size)
	R_break = np.average(med_Rs, weights = 1. / (np.array(gradients)))
	print R_break
	# plt.plot(med_Rs, gradients)
	# plt.show()
	return R_break, gradients, med_Rs

def define_boundaries(profile, window, chunk_size=2):
	R_break, gradients, med_Rs = find_break(profile, window, chunk_size)
	inner, outer = gradients[:R_break], gradients[R_break:]
	means = map(np.mean, [inner, outer])
	stds = map(np.std, [inner, outer])





if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 180
	target, info = tables[N], header.loc[N]
	result = F.fit_bulge_disc(target, info)
	lm.report_fit(result.params, show_correl=False)
	pos, delta = find_bulge_cutoff(target, info, result)

	brk, grads, grad_Rs = find_break(target, [pos, -info.sky_pos])

	fig, ax, ax2 = F.plot_basic(result, target, info)
	ax.axvline(x=target.R.values[pos])
	ax.axvline(x=target.R.values[-info.sky_pos])
	ax2.axvline(x=target.R.values[pos])
	ax2.axvline(x=target.R.values[-info.sky_pos])
	ax2.axvline(x=brk)
	ax.axvline(x=brk)
	plt.show()
