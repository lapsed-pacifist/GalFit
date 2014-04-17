import pandas as pd
import storage as S
import matplotlib.pyplot as plt
DIR = 'repository'


if __name__ == '__main__':
	T = pd.read_table(DIR+'\\fullcas_sky.dat', skiprows=0,skipinitialspace=True, escapechar='#',delimiter='\t')
	store = pd.HDFStore('store_again_100.h5', 'r')
	tables, info = S.import_directory()
	fits, info = store.fits, store.info
	tot = pd.merge(info, fits)
	tot = pd.merge(tot, T, on='ID')
	mask = ~tot.scale.isnull()
	plt.plot(tot.r_clust[mask], tot.BD_ratio[mask])
	plt.show()