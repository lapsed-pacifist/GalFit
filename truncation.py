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
	if len(diff) < 4:
		return 0, 0, False
	return diff[0], bulge[diff[0]] - disc[diff[0]], True

def find_sky_cutoff(profile, infoDF, sigma=3.):
	sky_cutoff = sigma * infoDF.sky_unc
	print sky_cutoff

def clean_profile(profile, window, spline_smooth=0, median_smooth_perc=0.05):
	sl = slice(window[0], window[1])
	mask = ~profile.M[sl].isnull().values
	y = profile.M[sl][mask].values
	x = profile.R[sl][mask].values
	tck = interpolate.splrep(x,y,s=spline_smooth)
	xnew = np.linspace(x[0],x[-1], len(x))
	f_size = len(x) * median_smooth_perc
	if f_size < 2:
		f_size = 2
	ynew = filt(interpolate.splev(xnew,tck,der=0), f_size)
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
	# plt.plot(med_Rs, gradients)
	# plt.show()
	if means[0] >= means[1]: #is inner > outer?
		inner_bound = np.where(inner[::-1] > means[0] + stds[0])
		outer_bound = np.where(outer < means[1] - stds[1])
	else:
		inner_bound = np.where(inner[::-1] < means[0] - stds[0])
		outer_bound = np.where(outer > means[1] + stds[1])
	try:
		inner_bound = med_Rs[brk_pos-1::-1][inner_bound[0][0]]
	except IndexError:
		inner_bound = med_Rs[0]
	try:
		outer_bound = med_Rs[brk_pos:][outer_bound[0][0]]
	except IndexError:
		outer_bound = med_Rs[-1]

	# plt.plot(med_Rs, gradients)
	# plt.axvline(R_break, linestyle='-')
	# colours = ['r', 'g']
	# for i in range(2):
	# 	plt.axhline(y=means[i]+stds[i], linestyle='--', color=colours[i])
	# 	plt.axhline(y=means[i]-stds[i], linestyle='--', color=colours[i])
	# 	plt.axhline(y=means[i], linestyle='-', color=colours[i])
	# plt.show()
	R_window = [profile.R.values[w] for w in window]
	return [R_window[0], inner_bound, outer_bound, R_window[1]], R_break, gradients, med_Rs

def fit_truncated(profile, infoDF, break_bounds, break_R):
	straight = lambda p, x,y,w: (y - (p['mu'].value +  (1.086 * x / p['h'].value))) / w
	R, M, M_err_up = profile.R.values, profile.M.values, profile.M_err_up.values
	mask_inner = (R <= break_R) & (R >= break_bounds[0]) & (~np.isnan(M_err_up))
	mask_outer = (R >= break_R) & (R <= break_bounds[-1]) & (~np.isnan(M_err_up))
	Rs = [R[mask_inner], R[mask_outer]]
	Ms = [M[mask_inner], M[mask_outer]]
	M_errs = [M_err_up[mask_inner], M_err_up[mask_outer]]
	P = lm.Parameters()
	P.add_many(('mu', 1., True, 10., 30.), ('h', 1., True, 0.05, 30.))
	ret = []
	for i in range(2):
		try:
			result = lm.minimize(straight, P, args=(Rs[i], Ms[i], M_errs[i]))
			result = np.array([result.params['mu'].value, result.params['h'].value])
		except TypeError:
			result = np.array([np.nan, np.nan])
		ret.append(result)
	return ret


def plot_truncation(axis, fit_results, overlap_perc=0.05):
	lims = axis.get_xlim()
	straight = lambda p, x: p[0] +  (1.086 * x / p[1])
	overlap = (lims[1] - lims[0]) * 0.5 * overlap_perc
	# print (lims[1] - lims[0]) * 0.5 * overlap_perc
	get_overlapx = lambda r, h: np.sqrt(r * r * h * h / ((h * h) + (1.086**2.)))
	overlap_in_x = [get_overlapx(overlap, p[1]) for p in fit_results]
	# print overlap_in_x
	break_R = (fit_results[0][0] - fit_results[1][0]) / ((1.086/fit_results[1][1]) - (1.086/fit_results[0][1]))
	# print fit_results
	# print break_R
	pltR1 = np.linspace(lims[0], break_R+(overlap_in_x[0]), 1000)
	pltR2 = np.linspace(break_R-(overlap_in_x[1]), lims[1], 1000)
	axis.plot(pltR1, straight(fit_results[0], pltR1), 'b-', linewidth=2)
	axis.plot(pltR2, straight(fit_results[1], pltR2), 'b-', linewidth=2)

def trunc_model(fit_results, r_brk, R):
	straight = lambda p, x: p[0] +  (1.086 * x / p[1])
	lowR, highR = R[R<r_brk], R[R>=r_brk]
	return np.append(straight(fit_results[0], lowR), straight(fit_results[1], highR))

def get_chi_trunc(normal_params, fit_results, r_brk, Y, R, W, zp):
	D = trunc_model(fit_results, r_brk, R)
	B = F.convert_mag(F.sersic(normal_params, zp, R, comp=True)[0], zp)
	return F.redchi(Y, B+D, W, 1.)

def fit_truncation(profile, infoDF, fit_result):
	pos, delta, success = find_bulge_cutoff(profile, infoDF, fit_result)
	if not success:
		return np.nan, np.nan, np.nan, np.nan, np.nan, False
	bounds, brk, grads, grad_Rs = define_boundaries(profile, [pos, -1])
	fits = fit_truncated(profile, infoDF, bounds, brk)
	return bounds, brk, pos, delta, fits, True


if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 72
	target, info = tables[N], header.loc[N]
	result = F.fit_bulge_disc(target, info)
	lm.report_fit(result.params, show_correl=False)
	fit_bound, brk_R, brk_pos, deltaB, fit_pair, success = fit_truncation(target, info, result)


	# import bootstrapping as B
	# print B.bootstrap_trunc(fit_pair, target.R, target.M, target.M_err_up, brk_pos, size=10)



	fig, ax, ax2 = F.plot_basic(result, target, info)
	for i in [ax, ax2]:
		i.axvline(x=brk_R)
		for b in fit_bound:
			i.axvline(x=b, linestyle='--')
	plot_truncation(ax, fit_pair)
	ax.set_xlim(target.R.values[0], target.R.values[-1])
	plt.show()
