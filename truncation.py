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

def exp_Imod(mu0, h, x, zp):
	I0 = 10 ** ((zp - mu0) / 2.5)
	return I0 * np.exp(-1.* x / h)

def mu_mod(mu0, h, x, zp):
	return mu0 + (1.086 * x / h)

def trunc_Imod(P, x, zp):
	innerR = x[x<P['Rbr'].value]
	outerR = x[x>=P['Rbr'].value]
	innerI = exp_Imod(P['mu01'].value, 1./P['invh1'].value, innerR, zp)
	outerI = exp_Imod(P['mu02'].value, P['h2'].value, outerR, zp)
	return np.append(innerI, outerI)

def trunc_mod(P, x, zp):
	innerR = x[x<P['Rbr'].value]
	outerR = x[x>=P['Rbr'].value]
	innerM = mu_mod(P['mu01'].value, 1./P['invh1'].value, innerR, zp)
	outerM = mu_mod(P['mu02'].value, P['h2'].value, outerR, zp)
	return np.append(innerM, outerM)


def fit_truncated(profile, infoDF, break_bounds, break_R, fix_brk=False):
	R, I, W = profile.R.values, profile.M.values, profile.M_err_down.values
	mask = (~np.isnan(I)) & (R > break_bounds[0])
	R = R[mask]
	I = I[mask]
	W = W[mask]
	res = lambda P, x,y,w,z: (y - trunc_mod(P, x, z)) / w
	
	P = lm.Parameters()
	hmin, hmax = 0.01, 20.0
	mumin, mumax = 14.0, 30.0

	P.add('mu02', value=25., vary=True, min=mumin, max=mumax)
	P.add('h2', value= 5., min=hmin, max=hmax)
	P.add('Rbr', value = break_R, vary=not fix_brk, min=break_bounds[1], max=break_bounds[2])
	P.add('deltah', value=1./300, max=(1./hmin)-(1./hmax), min=(-1./hmin)+(1./hmax))
	P.add('deltamu', expr='1.086 * Rbr * deltah')
	P.add('invh1', expr='(1./h2) + deltah')
	P.add('mu01', expr='mu02 - deltamu')
	result = lm.minimize(res, P, args=(R, I, W, infoDF.zp))
	innerresult = [P['mu01'].value, 1./P['invh1'].value]
	outerresult = [P['mu02'].value, P['h2'].value]
	# lm.report_fit(result.params, show_correl=False)
	return np.array([innerresult, outerresult]), P['Rbr'].value

def plot_truncation(axis, fit_results, break_R):
	lims = axis.get_xlim()
	straight = lambda p, x: p[0] +  (1.086 * x / p[1])
	# overlap = (lims[1] - lims[0]) * 0.5 * overlap_perc
	# get_overlapx = lambda r, h: np.sqrt(r * r * h * h / ((h * h) + (1.086**2.)))
	# overlap_in_x = [get_overlapx(overlap, p[1]) for p in fit_results]
	# break_R = (fit_results[0][0] - fit_results[1][0]) / ((1.086/fit_results[1][1]) - (1.086/fit_results[0][1]))

	pltR1 = np.linspace(lims[0], break_R+10., 1000)
	pltR2 = np.linspace(break_R-10., lims[1], 1000)
	axis.plot(pltR1, straight(fit_results[0], pltR1), 'r-', linewidth=1)
	axis.plot(pltR2, straight(fit_results[1], pltR2), 'b-', linewidth=1)

def trunc_model(fit_results, r_brk, R):
	straight = lambda p, x: p[0] +  (1.086 * x / p[1])
	lowR, highR = R[R<r_brk], R[R>=r_brk]
	return np.append(straight(fit_results[0], lowR), straight(fit_results[1], highR))

def get_chi_trunc(normal_params, fit_results, r_brk, Y, R, W, zp):
	D = trunc_model(fit_results, r_brk, R)
	B = F.convert_mag(F.sersic(normal_params, zp, R, comp=True)[0], zp)
	return F.redchi(Y, B+D, W, 1.)

def fit_truncation(profile, infoDF, fit_result):
	try:
		pos, delta, success = find_bulge_cutoff(profile, infoDF, fit_result)
	except TypeError:
		return 0.0, 0.0, 0.0, 0.0, 0.0, False
	if not success:
		return 0.0, 0.0, 0.0, 0.0, 0.0, False
	bounds, brk, grads, grad_Rs = define_boundaries(profile, [pos, -1])
	fits, R_brk = fit_truncated(profile, infoDF, bounds, brk)
	return bounds, R_brk, pos, delta, fits, True


if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 128
	target, info = tables[N], header.loc[N]
	print info.ID
	# target.I = (target.i_cts - (info.sky_unc)) / info.scale / info.scale
	# target.M = info.zp - (2.5 * np.log10(target.I))
	# up = target.I + target.I_err
	# down = target.I - target.I_err
	# Mdown = info.zp - (2.5 * np.log10(abs(up)))
	# Mup = info.zp - (2.5 * np.log10(abs(down)))
	# target['M_err_down'] = abs(target.M - Mdown)
	# target['M_err_up'] = abs(Mup - target.M)

	result = F.fit_bulge_disc(target, info)
	lm.report_fit(result.params, show_correl=False)
	fit_bound, brk_R, brk_pos, deltaB, fit_pair, success = fit_truncation(target, info, result)
	print fit_pair


	# import bootstrapping as B
	# print B.bootstrap_trunc(fit_pair, target.R, target.M, target.M_err_up, brk_pos, size=10)



	fig, ax, ax2 = F.plot_basic(result, target, info)
	for i in [ax, ax2]:
		i.axvline(x=brk_R)
		for b in fit_bound:
			i.axvline(x=b, linestyle='--')
	plot_truncation(ax, fit_pair, brk_R)
	ax.set_xlim(target.R.values[0], target.R.values[-1])
	plt.show()


