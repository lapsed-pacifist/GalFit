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

def extract_params(p):
	d = {name: par.value for name, par in p.iteritems()}
	return d

def batch_fit(table_list, infoDF, HDF, boot_size=1000, i_start=0, load_bar=True):
	L = loadbar.loadbar(len(table_list))
	for index,t in enumerate(table_list):
		i = index + i_start
		L.time()
		normal_fit_result = F.fit_bulge_disc(t, infoDF.loc[i])
		pars = extract_params(normal_fit_result.params)
		pars.update({'flag':normal_fit_result.success})

		boundaries, brk_R, brk_pos, deltaB, trunc_pair, success = T.fit_truncation(t, infoDF.loc[i], normal_fit_result)
		if not success:
			t_vars = {'brk_R':np.nan, 'brk_pos':np.nan, 'deltaB':np.nan,'inner_M':np.nan, 'inner_Re':np.nan, 'outer_M':np.nan, 'outer_Re':np.nan}
			bound_d = {'b'+str(i):np.nan for i in range(4)}
			t_vars.update(bound_d)
		else:
			inner, outer = trunc_pair
			bound_d = {'b'+str(i):v for i,v in enumerate(boundaries)}
			t_vars = {'brk_R':brk_R, 'brk_pos':brk_pos, 'deltaB':deltaB,'inner_M':inner[0], 'inner_Re':inner[1], 'outer_M':outer[0], 'outer_Re':outer[1]}
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
		HDF['bootstraps'] = wp

# def batch_bootstrap(table_list, infoDF, HDF, size=1000):
# 	L = loadbar.loadbar(len(table_list))
# 	for i,t in enumerate(table_list):
# 		L.time()
# 		P = HDF['fits'].copy()
# 		DF = B.bootstrap(P, t, infoDF.loc[i], size, False)
# 		if i == 0:
# 			HDF['bootstrap'] = pd.Panel([DF])
# 		else:
# 			HDF['bootstrap'][i] = DF
# 		L.progress()



if __name__ == '__main__':
	tables, info = S.import_directory()
	store = pd.HDFStore('store_again_100.h5')
	if len(store.keys()) != 0:
		print "Warning, store already exists! Data will be overwritten!"
		raw_input("to continue press enter...")
	select = tables[:]
	batch_fit(select, info, store, boot_size=100, load_bar=True)
	store['info'] = info
	print store.truncations.describe()
