import pandas as pd
import bootstrapping as B

store = pd.HDFStore('store_again.h5', 'r')
fits = store.fits
info = store.info
boot = store.bootstraps

fits = pd.concat([fits, info[['ID','cam','ax']]], axis=1)
for i, b in boot.iteritems():
	stds, mean = B.find_stats(b) #series
	if i == 0:
		meanDF = pd.DataFrame(mean).T
		stdupDF = pd.DataFrame(stds[1]).T
		stddownDF = pd.DataFrame(stds[0]).T
	else:
		meanDF = meanDF.append(mean, ignore_index=True)
		stdupDF = stdupDF.append(stds[1], ignore_index=True)
		stddownDF = stddownDF.append(stds[0], ignore_index=True)

print meanDF.loc[0]
print stdupDF.loc[0]
print stddownDF.loc[0]