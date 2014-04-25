import pandas as pd
import matplotlib.pyplot as plt

def graph_binned(r1, r2, r3, propx, propy, ax):
	y = [r1.propx]
	ax.plot()

if __name__ == '__main__':
	store_trunc = pd.HDFStore('fixed_truncations.h5')
	total = store_trunc['mainDF']
	r_regions = [35., 68.]
	R1 = total[total.r_clust < r_regions[0]]
	R2 = total[(total.r_clust > r_regions[0]) & (total.r_clust < r_regions[1])]
	R3 = total[total.r_clust > r_regions[1]]

	# fig = plt.figure()
	# ax = fig.add_subplot(111