import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import truncation as T
import storage as S
import lmfit as lm
import fit as F
from scipy import stats
import loadbar
import sys
pd.set_option('io.hdf.default_format','table')

def total_model(p, x, zp, comp=False):
	mue, Re, n = p['MB'].value, p['ReB'].value, p['nB'].value
	b = F.get_b_n(n)
	Ie = F.convert_I(mue, zp)
	bulge = Ie * np.exp(-1.* b * ( ((x / Re) ** (1 / n)) - 1 ))

	innerR = x[x<p['Rbr'].value]
	outerR = x[x>=p['Rbr'].value]
	innerI = T.exp_Imod(p['mu01'].value, p['h1'].value, innerR, zp)
	outerI = T.exp_Imod(p['mu02'].value, p['h2'].value, outerR, zp)
	disk = np.append(innerI, outerI)

	if comp:
		return bulge, innerI, outerI, innerR, outerR
	else:
		return disk + bulge

def test_trunc(profile, truncDF, infoDF, sigma=2):
	"""raises and lowers the sky by sigma"""
	bounds = truncDF[['b0', 'b1', 'b2', 'b3']].values
	sky = infoDF.sky_unc * sigma

	profile.I = (profile.i_cts - sky) / infoDF.scale / infoDF.scale
	profile.M = infoDF.zp - (2.5 * np.log10(profile.I))
	up = profile.I + profile.I_err
	down = profile.I - profile.I_err
	Mdown = info.zp - (2.5 * np.log10(abs(up)))
	Mup = info.zp - (2.5 * np.log10(abs(down)))
	profile['M_err_down'] = abs(profile.M - Mdown)
	profile['M_err_up'] = abs(Mup - profile.M)
	parsdown, brk = T.fit_truncated(profile, infoDF, bounds, truncDF.brk_R, fix_brk=True)

	profile.I = (profile.i_cts + sky) / infoDF.scale / infoDF.scale
	profile.M = infoDF.zp - (2.5 * np.log10(profile.I))
	up = profile.I + profile.I_err
	down = profile.I - profile.I_err
	Mdown = info.zp - (2.5 * np.log10(abs(up)))
	Mup = info.zp - (2.5 * np.log10(abs(down)))
	profile['M_err_down'] = abs(profile.M - Mdown)
	profile['M_err_up'] = abs(Mup - profile.M)
	parsup, brk = T.fit_truncated(profile, infoDF, bounds, truncDF.brk_R, fix_brk=True)

	return parsup, parsdown, np.array([[truncDF.inner_M, truncDF.inner_h], [truncDF.outer_M, truncDF.outer_h]]), brk, bounds

def get_chiKS(profile, fitsDF, mu_list, h_list, break_R, bounds, infoDF):
	P = lm.Parameters()
	P.add_many(('mu01',mu_list[0]), ('mu02', mu_list[1]), ('h1', h_list[0]), ('h2', h_list[1]))
	P.add('Rbr', break_R)
	P.add_many(('nB', fitsDF.nB), ('ReB', fitsDF.ReB), ('MB', fitsDF.MB), ('deltaRe', 1., True), ('BD_ratio', fitsDF.BD_ratio), ('ReD', fitsDF.ReD))
	res = (profile.I - total_model(P, profile.R.values, infoDF.zp, False)) / profile.I_err
	res_norm = F.sersic(P, infoDF.zp, profile.R.values, profile.I.values, profile.I_err, False) 
	return np.sum(res**2.), stats.kstest(res, 'norm'), np.sum(res_norm**2.), stats.kstest(res_norm, 'norm')

def compare_wsky(profile, fitsDF, truncDF, infoDF, sigma=2):
	up, down, norm, brk, bounds = test_trunc(profile, truncDF, infoDF, sigma)
	trunc_chi, trunc_KS, class_chi, class_KS = get_chiKS(profile, fitsDF, norm[:,0], norm[:,1], brk, bounds, infoDF)
	print bounds
	print type(bounds)
	# return (trunc_KS[1] - class_KS[1]) / trunc_KS[1]
	return down[1,1], norm[1,1], fitsDF.ReD/1.678, up[1,1]

def classify(profile, fitsDF, truncDF, infoDF):
	up, down, norm, brk, bounds = test_trunc(profile, truncDF, infoDF, 1.)
	h_fits = (fitsDF.ReD)/1.678
	if (h_fits < down[1,1]): #upbend
		cl = 1
	elif (h_fits > up[1,1]):#downbend
		cl = -1
	else:
		cl = 0
	# if ((norm[0,1] - norm[1,1]) / norm[0,1]) > 
	up_un, down_un, norm_un = [(i[0,0], i[0,1], i[1,0], i[1,1]) for i in (up, down, norm)]
	return up_un, down_un, norm_un, brk, bounds, cl

def condense(profile_table, fitsDF, truncDF, infoDF):
	DA = np.zeros([len(profile_table), 19]) # Data array for parameters
	L = loadbar.loadbar(len(profile_table))
	fails = 0
	for i, profile in enumerate(profile_table):
		L.time()
		try:
			plus, minus, trunc, brk, bounds, cl = classify(profile, fitsDF.loc[i], truncDF.loc[i], infoDF.loc[i])
			DA[i,:13] = plus + minus + trunc + (brk,)
			DA[i, 13:-2] = bounds
			DA[i,-2] = cl
			DA[i,-1] = i
		except TypeError:
			DA[i] = [np.nan]*19
			fails += 1
		L.progress()
		sys.stdout.write('fails=%i     ' % (fails))
	base = ['M1', 'h1', 'M2', 'h2']
	base = [[i+n for i in base] for n in ['_up', '_down','']]
	base = base[0]+base[1]+base[2]
	col_names = base + ['brk'] + ['b'+str(i) for i in range(1,5)] + ['trunc_type', 'index']
	return pd.DataFrame(DA, columns=col_names)

def total_classify(total_truncDF):
	


if __name__ == '__main__':
	tables, header = S.import_directory()
	N = 72
	target, info = tables[N], header.loc[N]
	store = pd.HDFStore('store_normal.h5')
	normal = store['truncations']
	fits = store['fits']
	normal = pd.concat([normal, header[['ID','cam','ax']]], axis=1)
	# tot = condense(tables, fits, normal, header)
	# store['total_truncations'] = tot
	tot = store['total_truncations']
	tot = pd.concat([tot, header[['ID','cam','ax']]], axis=1)
	tot = tot.set_index(['ID', 'ax', 'cam'])
	print tot.ix[1237665427553583111]



	# a = tot.trunc_type.hist()
	# plt.show()