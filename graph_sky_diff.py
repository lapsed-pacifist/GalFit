import storage as S
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
rc('font',**{'family':'serif','serif':['Times New Roman'], 'size':17})

if __name__ == '__main__':
	tables, header = S.import_directory()
	# fig = plt.figure()
	# ax = fig.add_axes([0, 1])
	# ax.plot(header.sky_level-header.sky, 'b.')
	# ax.hist(header.sky_level-header.sky,40,range=(-20,20),orientation='horizontal')
	a = (header.sky_level-header.sky)/ header.sky
	b = header.sky_unc - np.sqrt(header.sky) /  np.sqrt(header.sky)
	a.hist()
	print a.describe()
	print b.describe()
	plt.show()