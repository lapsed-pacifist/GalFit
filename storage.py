import pandas as pd
import glob
import os
import numpy as np
from copy import deepcopy
import sky_detect as sky
pd.set_option('io.hdf.default_format','table')
REP = 'repository/'

def import_ascii(name):
	with open(name, 'r') as f:
		Head = f.readline()[5:].rstrip().split(' ')
		Head.append(os.path.basename(name).split('_')[0])
	T = pd.read_table(name, skiprows=1,skipinitialspace=True, escapechar='#',delimiter=' ')
	T = T.rename(columns=lambda x: x.replace(' ', ''))
	n = ['ID', 'scale', 'zp', 'sky', 'cam']
	dtypes = [int, float, float, float, str]
	Head = [t(x) for t,x in zip(dtypes,Head)]
	H = {n[i]: v for i,v in enumerate(Head)}
	basis = T.ix[:,'method':'r_pixel']
	axis1 = pd.concat([basis,T.ix[:,'n1':'i1_cts_err']], axis=1)
	axis2 = pd.concat([basis,T.ix[:,'n2':'i2_cts_err']], axis=1)
	axis1 = axis1.rename(columns=lambda x: x.replace('1', ''))
	axis2 = axis2.rename(columns=lambda x: x.replace('2', ''))
	H1, H2 = H.copy(), H
	H1['ax'], H2['ax'] = 1, 2
	return axis1, axis2, H1, H2 

def import_directory(direct=REP, exclusions='exclusions.txt'):
	with open(REP+'exclusions.txt', 'r') as f:
		ex = f.read().split('\n')
	files = [fn for fn in glob.glob(REP+'*.ascii')]
	filtered = [f for f in files if os.path.basename(f).split('_')[1] not in ex]
	DF_list = [import_ascii(i) for i in filtered]
	table_list, infoDF = construct_info(DF_list)

	sky_vals = map(list, zip(*[sky.detect_sky(table.i_cts.values) for table in table_list]))
	infoDF['sky_level'], infoDF['sky_unc'], infoDF['sky_pos'] = sky_vals
	for i, v in enumerate(table_list):
		v.i_cts = v.i_cts - infoDF.sky_level.loc[i]

	table_list = convert(table_list, infoDF)
	table_list = get_errors(table_list, infoDF)
	table_list = convert_I_err(table_list, infoDF)
	return table_list, infoDF


def construct_info(table_list):
	dict_list = zip(*table_list)
	return dict_list[0]+dict_list[1], pd.DataFrame(list(dict_list[2]+dict_list[3]))

def convert(table_list, infoDF):
	for i, T in enumerate(table_list):
		H = infoDF.loc[i]
		T['R'] = T.r_pixel * H['scale']
		T['I'] = (T.i_cts) / (H['scale']**2.)
		T['M'] = H['zp'] - (2.5 * np.log10(T.I))
	return table_list

def convert_I_err(table_list, infoDF):
	for i, T in enumerate(table_list):
		up = T.I + T.I_err
		down = T.I - T.I_err
		Mdown = infoDF.zp[i] - (2.5 * np.log10(abs(up)))
		Mup = infoDF.zp[i] - (2.5 * np.log10(abs(down)))
		T['M_err_down'] = abs(T.M - Mdown)
		T['M_err_up'] = abs(Mup - T.M)
	return table_list

def get_errors(table_list, infoDF):
	for i, T in enumerate(table_list):
		if infoDF.cam[i] == 'sdss':
			t, G = 54., 7.43
		elif infoDF.cam[i] == 'mega':
			t, G = 35., 1.7
		else: t, G = 1.,1.
		phot_err = np.sqrt(((T.i_cts + infoDF.sky) / 1.) * G)
		i_err = phot_err * 1. / (G)
		I_err = (i_err) / (infoDF.scale[i] ** 2.)
		total = np.sqrt((I_err **2.) + ((T.i_cts_err / (infoDF.scale[i]**2.)) **2.) + (infoDF.sky_unc[i] / (infoDF.scale[i]**2.)))
		# total = [i if i > 1000. else 1000. for i in total]
		T['I_err'] = total
	return table_list



if __name__ == '__main__':
	tables, info = import_directory()
	# for i, t in enumerate(tables):
	# 	print info.loc[i].ax, info.loc[i].cam
	

	# print tables[0][['I', 'I_err', 'M','M_err_down', 'M_err_up']].head(18)
	import matplotlib.pyplot as plt
	N = 20
	plt.errorbar(tables[N].R.values, tables[N].M.values, yerr=[tables[N].M_err_down.values, tables[N].M_err_up.values], fmt='b.')
	plt.axvline(x=tables[N].R.values[-info.sky_pos[N]], linestyle='--', color='r')
	# plt.axhline(y=0)
	plt.title(str(info.ID[N])+info.cam[N]+str(info.ax[N]))
	plt.ylim([35,15])
	plt.show()

	# plt.plot(tables[0].R, tables[0].M, 'b.')
	# print info.loc[0].sky_level
	# plt.ylim([35,15])
	# plt.show()