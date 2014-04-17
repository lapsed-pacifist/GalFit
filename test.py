import pandas as pd
pd.set_option('io.hdf.default_format','table')

store = pd.HDFStore('store_again_100.h5')
print store.bootstraps[74].describe()