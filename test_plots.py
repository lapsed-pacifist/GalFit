import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# f, axarr = plt.subplots(1, 2)
# f.subplots_adjust(wspace=0., hspace=0.)


# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.

# divider = make_axes_locatable(ax)
# ax2 = divider.append_axes("right", size="100%", pad=0.00)
# ax3 = divider.append_axes("bottom", size="100%", pad=0.00)
# divider2 = make_axes_locatable(ax3)
# ax4 = divider2.append_axes("right", size="100%", pad=0.00)
# cax.set_yticklabels(cax.get_yticklabels(), visible=False)


# plt.colorbar(im, cax=cax)
# plt.show()
f = plt.figure()
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator


