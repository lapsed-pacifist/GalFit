import time
import storage as S
import numpy as np
import matplotlib.pyplot as plt

def convert_mag(I, zp):
	return zp - (2.5 * np.log10(I))

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

def find_mu_crit(profile, infoDF, delta=0.2):
	"""gets the critical SB for a corrected profile, where the +-sigmas
	differ by delta"""
	plus, minus = profile.I.values + infoDF.sky_unc, profile.I.values - infoDF.sky_unc
	diff = np.abs(convert_mag(minus, infoDF.zp) - convert_mag(plus, infoDF.zp))
	return profile.M.values[diff >  delta][0]



if __name__ == '__main__':
	tables, header = S.import_directory()
	i=2
	target, info = tables[i], header.loc[i]
	print len(tables)
	A = header.sky - header.sky_level
	print A.describe()
	crit = find_mu_crit(target, info)

	# select = [10]
	# for i in select:
	# 	target, info = tables[i], header.loc[i]
	# 	mean, std, cutoff = detect_sky(target.i_cts.values)
	# 	# print mean, std, cutoff 

	plt.plot(target.R.values, target.M.values, 'b.')
	plt.axhline(y=crit)
	# plt.axhline(y=0)
	# plt.axhline(y=0-std)
	# plt.axhline(y=0+std)
	# plt.axvline(x=target.R.values[-cutoff])
	# plt.title(str(info.ID)+info.cam+str(info.ax))
	plt.show()