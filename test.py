import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import storage as S

def hist(x, n_bin, axis):
	hist, bins = np.histogram(x, bins=50)
	width = 1. * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	axis.bar(center, hist, align='center', width=width)

if __name__ == '__main__':	

	store = pd.HDFStore('store_trunc_boot200.h5', 'r')
	boot = store.trunc_boot
	print boot[0].