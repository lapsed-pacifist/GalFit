import fit as F
import storage as S
import bootstrapping as B
import truncation as T
import sky_detect as sky
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loadbar
pd.set_option('io.hdf.default_format','table')
import sys
import winsound
import lmfit as lm

def extract_params(p):
	d = {name: par.value for name, par in p.iteritems()}
	return d

def batch_fit(table_list, infoDF, HDF, boot_size=1000, i_start=0, load_bar=True):
	L = loadbar.loadbar(len(table_list))
	trunc_fail = 0
	for index,t in enumerate(table_list):
		i = index + i_start
		L.time()

		normal_fit_result = F.fit_bulge_disc(t, infoDF.loc[i])
		pars = extract_params(normal_fit_result.params)
		pars.update({'flag':normal_fit_result.success})
		try:
			boundaries, brk_R, brk_pos, deltaB, trunc_pair, success = T.fit_truncation(t, infoDF.loc[i], normal_fit_result)
		except TypeError:
			boundaries, brk_R, brk_pos, deltaB, trunc_pair, success = 0.0, 0.0, 0.0, 0.0, 0.0, False

		if not success:
			trunc_fail += 1

		if not success:
			t_vars = {'brk_R':np.nan, 'brk_pos':np.nan, 'deltaB':np.nan,'inner_M':np.nan, 'inner_h':np.nan, 'outer_M':np.nan, 'outer_h':np.nan, 'success':False}
			bound_d = {'b'+str(i):np.nan for i in range(4)}
			t_vars.update(bound_d)
		else:
			inner, outer = trunc_pair
			bound_d = {'b'+str(i):v for i,v in enumerate(boundaries)}
			t_vars = {'brk_R':brk_R, 'brk_pos':brk_pos, 'deltaB':deltaB,'inner_M':inner[0], 'inner_h':inner[1], 'outer_M':outer[0], 'outer_h':outer[1], 'success':False}
			t_vars.update(bound_d)

		DF_boot = B.bootstrap(normal_fit_result.params, t, infoDF.loc[i], boot_size, False)

		if i == 0:
			HDF['fits'] = pd.DataFrame(pars, index=[0], dtype=float)
			HDF['truncations'] = pd.DataFrame(t_vars, index=[0], dtype=float)
			# HDF['bootstraps'] = pd.Panel({0:DF_boot})
			wp = pd.Panel({0:DF_boot})
		else:
			HDF['fits'] = HDF['fits'].append(pars, ignore_index=True)
			HDF['truncations'] = HDF['truncations'].append(t_vars, ignore_index=True)
			# HDF['bootstraps'][i] = DF_boot
			wp[i] = DF_boot
		if load_bar:
			L.progress()
			sys.stdout.write('fails=%i   ' % trunc_fail)
		HDF['bootstraps'] = wp

def batch_truncation(table_list, infoDF, truncDF, fitDF, HDF, boot_size=200, i_start=0, load_bar=True):
	try:
		if load_bar : L = loadbar.loadbar(len(table_list) * boot_size, boot_size)
		else: L = None
		for i,t in enumerate(table_list):
			index = i + i_start
			f = fitDF.loc[index]
			fP = lm.Parameters()
			fP.add_many(('MB', f.MB), ('ReB', f.ReB), ('deltaRe', 1., True),('nB', f.nB), ('BD_ratio', f.BD_ratio), ('ReD', f.ReD))
			DF_boot = B.bootstrap_sky(t, infoDF.loc[index], truncDF.loc[index], fP, boot_size, L)
			if index == 0:
				wp = pd.Panel({0:DF_boot})
			else:
				wp[index] = DF_boot
		HDF['trunc_boot'] = wp
	except:
		winsound.PlaySound('SystemHand', winsound.SND_ALIAS)


if __name__ == '__main__':
	tables, info = S.import_directory()
	store_name = 'fixed_truncations.h5'
	store = pd.HDFStore(store_name)
	truncs = store.truncations
	fits = store.fits
	# if len(store.keys()) != 0:
	# 	print "Warning, store %s already exists! Data will be overwritten!" % (store_name)
	# 	raw_input("to continue press enter...")
	# select = tables[:]
	# batch_fit(select, info, store, boot_size=1, load_bar=True)
	# store['info'] = info
	# print store.truncations.describe()
	
	batch_truncation(tables[:], info, truncs, fits, HDF=store, i_start=0, load_bar=True)
	print store['trunc_boot']

# print truncs.rename(columns={'inner_M':'mu_inner', 'outer_M':'mu_outer', 'inner_Re':'h_inner', 'outer_Re':'h_outer'})