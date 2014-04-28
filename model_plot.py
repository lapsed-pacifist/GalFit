import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman']})

def get_b_n(m):
	#   find sersic b_n coefficient'favoritecolor'
	#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
	if m == 0.0: return -0.3271
	b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
				+ 131.0/m/m/m/1148175.0  \
				- 2194697.0/m/m/m/m/30690717750.0
	return b_n

def sersic(Me, Re, n, X):
	A = ((X / Re) ** (1./n)) - 1
	return Me + (2.5 * get_b_n(n) * A / np.log(10))

def convert_I(mu, zp=30.):
	if type(mu) == list:
		mu = np.array(mu)
	return 10 ** ((zp - mu) / 2.5)

def convert_mag(I, zp=30.):
	return zp - (2.5 * np.log10(I))

def trunc_disc(M_inner, hs, brk, X):
	M_outer = (1.086 * brk * ((1./hs[0]) - (1./hs[1]))) + M_inner
	Ms = [M_inner, M_outer]
	disc = lambda mu, h, x:  float(mu) + (1.086 * x / h)
	Rs = [X[X < brk], X[X >=brk]]
	discs = [disc(Ms[i], hs[i], Rs[i]) for i in range(2)]
	return discs, Rs, disc(Ms[1], hs[1], brk)

def plot_profile(Me, Re, n, M_inner, hs, brk, sky, X, ax, title, box=True, text=True, leg=True):
	bulge = sersic(Me, Re, n, X)
	discs, disc_R, mu_brk = trunc_disc(M_inner, hs, brk, X)
	tot = convert_I(bulge) + convert_I(np.hstack((discs[0],discs[1]))) + sky

	ax.plot(X, convert_mag(tot), 'k--', label='total', linewidth=2)
	ax.plot(X, bulge, 'g:', label='bulge', linewidth=2)
	ax.plot(disc_R[0], discs[0], color='b', linestyle=':',  label='inner disc', linewidth=2)
	ax.plot(disc_R[1], discs[1], color='b', linestyle=':',  label='outer disc', linewidth=2)
	

	ax.set_ylim([35, 15])
	ax.set_ylabel('$\mu (R)$ [mag]', fontsize=14)
	ax.set_xlabel('R [arcsec]', fontsize=14)
	if text:
		titles = ['\\mu_e', 'R_e', 'n', '\\mu_{0, inner}', 'h_{inner}', 'h_{outer}', 'R_{brk}']
		var = [Me, Re, n, M_inner, hs[0], hs[1], brk]
	else:
		titles = ['h_{inner}', 'h_{outer}']
		var = hs
	s = ''
	if text is not None:
		for name, val in zip(titles, var):
			s += '$%s = %.1f$\n' % (name, val)
		ax.text(0.98, 0.98, s, transform=ax.transAxes, va='top', ha='right', fontsize=16)
	if box:
		ax.plot(brk, mu_brk, 'rs', markersize=10)
	if leg:
		ax.legend(loc=3, fancybox=True)
	ax.set_title(title)


if __name__ == '__main__':
	R = np.linspace(0,50,1000)
	fig = plt.figure()
	fig.set_facecolor('white')
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)
	plot_profile(20., 4., 0.7, 18., [5., 5.], 30., 0., R, ax1, 'Correct Sky', box=False, text=None, leg=False)
	plot_profile(20., 4., 0.7, 18., [5., 5.], 30., -10., R, ax2, 'Underestimated Sky', box=False, text=None, leg=False)
	plot_profile(20., 4., 0.7, 18., [5., 5.], 30., 10., R, ax3, 'Overestimated Sky', box=False, text=None, leg=False)
	# plot_profile(20., 4., 0.7, 18., [5., 3.], 30., 0., R, ax2, 'Type-II (downbending)', box=True, text=False, leg=False)
	# plot_profile(20., 4., 0.7, 18., [3., 5.], 30., 0., R, ax3, 'Type-III (upbending)', box=True, text=False, leg=False)
	
	
	plt.show()

