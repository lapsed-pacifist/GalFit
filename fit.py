import storage as S
from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt
import bisect
import lmfit as lm
from scipy.special import gamma

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
		front = nB * gamma(2*nB) * np.exp(bB) / (bB ** (2*nB))
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
			   ('deltaRe', 5., True, 0., 1000.),
			   ('BD_ratio', 0.5, False, 0.001, 1.))
	fitter = lm.Minimizer(sersic, P, fcn_args=(infoDF.zp, profile.R.values, profile.I.values, profile.I_err.values))
	fitter.leastsq()
	
	P['deltaRe'].value = 5.
	P['BD_ratio'].value = 0.5

	fix_params(P, ['BD_ratio', 'deltaRe'], True)
	fitter = lm.Minimizer(sersic, P, fcn_args=(info.zp, profile.R.values, profile.I.values, profile.I_err.values))
	fitter.leastsq()
	return fitter


if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 211
	target, info = tables[N], header.loc[N]
	import time
	s = time.clock()
	result = fit_bulge_disc(target, info)
	print time.clock() - s
	lm.report_fit(result.params, show_correl=False)


	fig = plt.figure()
	ax = fig.add_subplot(211)
	ax.errorbar(target.R, target.M, yerr=(target.M_err_down.values, target.M_err_up.values), fmt='b.')
	pltR = np.linspace(0,target.R.values[-1], 1000)
	bulge, disc = sersic(result.params, info.zp, pltR, comp=True, show=True)
	
	ax.plot(pltR, convert_mag(bulge, info.zp), 'g:')
	ax.plot(pltR, convert_mag(disc, info.zp), 'r--')
	ax.plot(pltR, convert_mag(bulge+disc, info.zp), 'k--')
	ax.set_title(str(info.ID) + str(info.cam) + str(info.ax))
	ax.set_ylim(35,15)

	ax2 = fig.add_subplot(212)
	bulge, disc = sersic(result.params, info.zp, target.R.values, comp=True, show=True)
	ax2.plot(target.R, target.M - convert_mag(bulge + disc, info.zp), 'b.')
	string = r'$\chi^{2}_{\nu} = $%.2f' % result.redchi
	ax2.text(0.98, 0.98, string, transform=ax2.transAxes, ha='right', va='top')
	plt.show()