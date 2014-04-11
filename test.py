import lmfit as lm
import numpy as np

def func(p, x, y):
	return (p['m'].value * x) - y

def my_func(x):
    return (x /3.)

P = lm.Parameters()
P.add('m', 10.)
P.add('L', value=None, vary=False, expr='my_func(1)')
X = np.arange(100)
Y = X*2.

# fitter = lm.Minimizer(func, P, fcn_args=(X, Y))
# fitter.asteval.symtable['my_func'] = my_func
# result = fitter.leastsq()
# lm.report_fit(fitter.params)

print 2. /2 /2