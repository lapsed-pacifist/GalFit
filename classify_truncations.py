import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import storage as S
import bootstrapping as B

def get_stats(bootP):
	for i,b in bootP.iteritems():
		stds, mean = B.find_stats(b)
		if i == 0:
			meanDF = pd.DataFrame(mean).T
			stdDFs = [pd.DataFrame(stds[0], index=[0]).T, pd.DataFrame(stds[0], index=[0]).T]
		else:
			meanDF = meanDF.append(mean, ignore_index=True)
			stdDFs[0] = stdDFs[0].append(stds[0], ignore_index=True)
			stdDFs[1] = stdDFs[1].append(stds[1], ignore_index=True)
	return stdDFs, meanDF

def streamline_DF(truncDF, bootP, infoDF):
	truncDF.rename(columns={'inner_M':'M1', 'inner_h':'h1', 'outer_M':'M2', 'outer_h':'h2'}, inplace=True)
	stdDFs, meanDF = get_stats(bootP)
	truncDF = pd.concat([truncDF, infoDF[['ID','cam','ax']]], axis=1)
	meanDF = pd.concat([meanDF, infoDF[['ID','cam','ax']]], axis=1)
	stdDFs[0] = pd.concat([stdDFs[0], infoDF[['ID','cam','ax']]], axis=1)
	stdDFs[1] = pd.concat([stdDFs[1], infoDF[['ID','cam','ax']]], axis=1)
	truncDF['accepted'] = True # accepted column
	for i, ID in info.ID.iteritems():
		row = truncDF[truncDF.ID == ID]
		t = row[['M2', 'h2']]
		sup = stdDFs[1][stdDFs[1].ID == ID][['M2', 'h2']]
		sdo = stdDFs[0][stdDFs[0].ID == ID][['M2', 'h2']]
		w = 2. / (sup + sdo)
		av = (w * t).sum() / w.sum()
		print 'average', av.h2
		print t.h2
		# for j, prof in t.iterrows():
			# if average within parameter +/- std
			# overlap = ((prof - sdo.loc[j]) <= (av + av_err)) & ((prof + sup.loc[j]) >= (av - av_err))
			# print overlap.M2 & overlap.h2

			# if ((prof.M2-(5. * sdo.loc[j].M2)) <= av.M2 <= (prof.M2 + (5. * sup.loc[j].M2)))\
			# 	and ((prof.h2-(5. * sdo.loc[j].h2)) <= av.h2 <= (prof.h2 + (5.* sup.loc[j].h2))) and (row.loc[j].isnull().sum() == 0):
			# else:
			# 	truncDF.accepted[j] = False
	return truncDF


if __name__ == '__main__':
	tables, info = S.import_directory()
	store = pd.HDFStore('store_normal.h5')
	truncations = store['truncations']
	fits = store['fits']
	boots = pd.HDFStore('store_trunc_boot.h5').trunc_boot
	print streamline_DF(truncations, boots, info)



