import time
import storage as S
import numpy as np
import matplotlib.pyplot as plt
from fit import convert_mag, sersic

def chunks(l, n, lower_bound=0):
    chnks = [l[i:i+n] for i in range(0, len(l), n)]
    if len(chnks[-1]) < lower_bound:
    	steal = lower_bound - len(chnks[-1])
    	chnks[-1] = np.append(chnks[-2][n-steal:], chnks[-1])
    return chnks

def get_stats(List, limit=0):
	mean, std = np.mean(List), np.std(List)
	if std < limit: std = float(limit)
	return mean, std

def find_change(mean_list, std_list):
	End = len(mean_list) - 1
	for i in range(len(mean_list)):
		if i == End:
			raise IndexError("End Reached, change not found")
		else:
			if mean_list[i+1] > mean_list[i] + (1*std_list[i]):
				if mean_list[i+2] > mean_list[i+1] + (1*std_list[i+1]):
					return i+1
				else:
					if i+1 == End:
						raise IndexError("End Reached, change not found")
					if mean_list[i+2] > mean_list[i] + (1*std_list[i]):
						return i+2
					else:
						continue
			else:
				continue

def find_cutoff(List, chunk_size, start_window_size=4, limit=0):
	main, start = List[start_window_size:], List[0:start_window_size]
	start_mean, start_std = get_stats(start, limit)
	means, stds = zip(*map(get_stats, chunks(main, chunk_size)))
	means, stds = np.insert(means, 0, start_mean), np.insert(stds, 0, start_std)
	boundary = find_change(means, stds)
	return (boundary*chunk_size) + start_window_size

def clip(List, sig=3, show=False):
	if type(List) is not np.ndarray:
		List = np.array(List)
	m, s = np.mean(List), np.std(List)
	return List[(List < m+(sig*s)) & (List > m-(sig*s))]

def clipped_stats(List, sig=3):
	X = clip(List, sig, show=True)
	return np.mean(X), np.std(X)

def detect_sky(I):
	pos = find_cutoff(I[::-1], 3, 3, 0)
	clipped = clip(I[-pos:])
	return np.mean(clipped), np.std(clipped), pos

def find_mu_crit(I_corrected, sky_error, zp, fit_params, delta=0.2):
	"""gets the critical SB for a corrected profile, where the +-sigmas
	differ by delta"""
	plus, minus = I_corrected + sky_error, I_corrected - sky_error
	diff = convert_mag(plus, zp) - convert_mag(minus, zp)




if __name__ == '__main__':
	tables, header = S.import_directory()
	print len(tables)

	select = [0]
	for i in select:
		target, info = tables[i], header.loc[i]
		mean, std, cutoff = detect_sky(target.i_cts.values)
		print mean, std, cutoff 

		lim = 25
		plt.plot(target.R.values, target.i_cts.values[:]-mean, 'b.')
		plt.axhline(y=0)
		plt.axhline(y=0-std)
		plt.axhline(y=0+std)
		plt.axvline(x=target.R.values[-cutoff])
		plt.title(str(info.ID)+info.cam+str(info.ax))
		plt.show()