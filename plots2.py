import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
rc('font',**{'family':'serif','serif':['Times New Roman']})

params={'axes.labelsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'figure.figsize':(8,6), 'figure.facecolor':'white'}
plt.rcParams.update(params)

def make_gridspec(f, MR_RATIO, left, right, top, bottom, hspace, hist=False, norm_res = False):
	gs = gridspec.GridSpec(MR_RATIO+1+norm_res, 6)
	gs.update(left=left, right=right, top=top, bottom=bottom, hspace=hspace)
	if norm_res:
		if hist:
			m = f.add_subplot(gs[0:MR_RATIO, 1:6])
			r = f.add_subplot(gs[MR_RATIO:MR_RATIO+1, 1:6], sharex=m)
			nr = f.add_subplot(gs[MR_RATIO+1:, 1:6], sharex=m)
		else:
			m = f.add_subplot(gs[0:MR_RATIO-2, :])
			r = f.add_subplot(gs[MR_RATIO-2:MR_RATIO-1, :], sharex=m)
			nr = f.add_subplot(gs[MR_RATIO-1:, :], sharex=m)
	else:
		m = f.add_subplot(gs[0:MR_RATIO-1, :])
		r = f.add_subplot(gs[MR_RATIO-1:, :], sharex=m)


def make_axes():
	fig = plt.figure()
	make_gridspec(fig, 6,0.05,0.98,0.98,0.05,0,True, True)
	plt.show()

if __name__ == '__main__':
	make_axes()