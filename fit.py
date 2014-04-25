import storage as S
import numpy as np
import matplotlib.pyplot as plt
import bisect
import lmfit as lm
from scipy.special import gamma
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman'], 'size':18})

def get_b_n(m):
	#   find sersic b_n coefficient'favoritecolor'
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def convert_mag(I, zp):
	return zp - (2.5 * np.log10(I))

def convert_I(M, zp):
	return 10. ** ((zp - M) / 2.5)

def sersic(p, zp, x, y=None, w=None, comp=False, show=False):
	"""use mue, Re and n"""
	if not p['deltaRe'].vary:
		mue, Re, n = p['MB'].value, p['ReB'].value, p['nB'].value
		b = get_b_n(n)
		Ie = convert_I(mue, zp)
		SB = Ie * np.exp(-1.* b * ( ((x / Re) ** (1 / n)) - 1 ))
	else:
		meB, ReB, nB, BD_ratio, ReD = p['MB'].value, p['ReB'].value, p['nB'].value, p['BD_ratio'].value, p['ReD'].value
		bB = get_b_n(nB)
		bD = get_b_n(1.)
		IeB = convert_I(meB, zp)
		try:
			front = nB * gamma(2*nB) * np.exp(bB) / (bB ** (2*nB))
		except ValueError:
			front = nB * gamma(2*nB) * np.exp(bB) / (-1.*(abs(bB) ** (2*nB)))
		h = ReD / bD
		IeD = front * ReB * ReB * IeB / BD_ratio / h / h / np.exp(bD)
		B = IeB * np.exp(-1.* bB * ( ((x / ReB) ** (1 / nB)) - 1 ))
		D = IeD * np.exp(-1.* bD * ( (x / ReD) - 1 ))
		SB = B + D
	if y is None:
		if comp == False:
			return SB
		else:
			return B, D
	else:
		if w is None: w = np.ones(y.shape)
		return (y - SB) / w

def fix_params(p, fixed, unfix=False):
	if type(fixed) is list:
		for n in fixed:
			p[n].vary = unfix
	else:
		for n, v in fixed.iteritems():
			p[n].vary = unfix
			if v is not None:
				p[n].value = v

def BDratio(MB, ReB, nB, MD, ReD):
	h = ReD / 1.678
	I_ratio = (10 ** ((MD - MB) / 2.5)) * np.exp(-1.*get_b_n(1.))
	scale_ratio = ReB * ReB / h / h
	front = BD_y(nB)
	return front * scale_ratio * I_ratio

def BD_y(m):
	b = get_b_n(m)
	return m * gamma(2*m) * np.exp(b) / (b ** (2*m))

def fit_bulge_disc(profile, infoDF):
	P = lm.Parameters()
	P.add_many(('MB', 20., True, 10., 30.),
			   ('ReB', 5., True, 0.1, 20.),
			   ('nB', 4., True, 0.01, 10.),
			   ('ReD', None, False, None, None, 'ReB + deltaRe'),
			   ('deltaRe', 5., True, 1., 20.),
			   ('BD_ratio', 0.5, False, 0.001, 1.))
	fitter = lm.Minimizer(sersic, P, fcn_args=(infoDF.zp, profile.R.values, profile.I.values, profile.I_err.values))
	fitter.leastsq()
	
	P['deltaRe'].value = 5.
	P['BD_ratio'].value = 0.5

	fix_params(P, ['BD_ratio', 'deltaRe'], True)
	fitter = lm.Minimizer(sersic, P, fcn_args=(infoDF.zp, profile.R.values, profile.I.values, profile.I_err.values))
	fitter.leastsq()
	# sb = sersic(P, infoDF.zp, profile.R.values, profile.I.values, profile.I_err.values, comp=False)
	# A = ((profile.I.values - sb) / profile.I_err.values)
	# print np.sum(A**2.) / (5+len(profile.I.values)), fitter.redchi
	# from scipy import stats
	# print stats.kstest(A, 'norm')
	
	return fitter


def redchi(y, model, weights, free_params):
	resid = (y - model / weights) **2.
	return np.sum(resid) / (len(y) + free_params)

def plot_basic(fit_result, profile, infoDF, ax=None):
	if ax is None:
		fig = plt.figure()
		fig.set_facecolor('white')
		ax = fig.add_subplot(111)
	ax.errorbar(profile.R, profile.M, yerr=(profile.M_err_down.values, profile.M_err_up.values), fmt='b.')
	pltR = np.linspace(0,profile.R.values[-1], 1000)
	bulge, disc = sersic(fit_result.params, infoDF.zp, pltR, comp=True, show=True)
	
	ax.plot(pltR, convert_mag(bulge, infoDF.zp), 'g:')
	ax.plot(pltR, convert_mag(disc, infoDF.zp), 'r--')
	ax.plot(pltR, convert_mag(bulge+disc, infoDF.zp), 'k--')
	ax.set_title(str(infoDF.ID) + str(infoDF.cam) + str(infoDF.ax))
	ax.set_ylim(35,15)

	# bulge, disc = sersic(fit_result.params, infoDF.zp, profile.R.values, comp=True, show=True)
	# ax2.plot(profile.R, (profile.M - convert_mag(bulge + disc, infoDF.zp)) / 1., 'b.')
	# string = r'$\chi^{2}_{\nu} = $%.2f' % fit_result.redchi
	# ax2.text(0.98, 0.98, string, transform=ax2.transAxes, ha='right', va='top')
	return ax#, ax2


if __name__ == '__main__':
	tables, header = S.import_directory()
	N_list = [10]
	for N in N_list:
		target, info = tables[N], header.loc[N] 
		result = fit_bulge_disc(target, info)
		lm.report_fit(result.params, show_correl=False)
	plot_basic(result, target, info)
	plt.show()
