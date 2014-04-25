import pandas as pd
import numpy as np
from matplotlib.pyplot import show

def find_stats(DF):
	mean = DF.mean(axis=0)
	upper = DF[DF >= mean]
	lower = DF[DF <= mean]
	return [sided_std(lower, mean), sided_std(upper, mean)], mean

def sided_std(df, mean):
	dev2 = (df - mean)**2.
	N = (1./(len(df))) 
	standard  = np.sqrt(dev2.sum(axis=0) * N)
	return standard

if __name__ == '__main__':
	a = pd.Series(range(5), index=list('abcde'))
	a.rename({'a':'aaaa'}, inplace=True)
