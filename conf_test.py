from scipy.stats import norm
from numpy import sqrt
def binP(N, p, x1, x2):
	p = float(p)
	q = p/(1-p)
	k = 0.0
	v = 1.0
	s = 0.0
	tot = 0.0

	while(k<=N):
			tot += v
			if(k >= x1 and k <= x2):
					s += v
			if(tot > 10**30):
					s = s/10**30
					tot = tot/10**30
					v = v/10**30
			k += 1
			v = v*q*(N+1-k)/k
	return s/tot

def calcBin(vx, vN, vCL = 95):
	'''
	Calculate the exact confidence interval for a binomial proportion

	Usage:
	>>> calcBin(13,100)    
	(0.07107391357421874, 0.21204372406005856)
	>>> calcBin(4,7)   
	(0.18405151367187494, 0.9010086059570312)
	''' 
	vx = float(vx)
	vN = float(vN)
	#Set the confidence bounds
	vTU = (100 - float(vCL))/2
	vTL = vTU

	vP = vx/vN
	if(vx==0):
			dl = 0.0
	else:
			v = vP/2
			vsL = 0
			vsH = vP
			p = vTL/100

			while((vsH-vsL) > 10**-5):
					if(binP(vN, v, vx, vN) > p):
							vsH = v
							v = (vsL+v)/2
					else:
							vsL = v
							v = (v+vsH)/2
			dl = v

	if(vx==vN):
			ul = 1.0
	else:
			v = (1+vP)/2
			vsL =vP
			vsH = 1
			p = vTU/100
			while((vsH-vsL) > 10**-5):
					if(binP(vN, v, 0, vx) < p):
							vsH = v
							v = (vsL+v)/2
					else:
							vsL = v
							v = (v+vsH)/2
			ul = v
	return (dl, ul)

def wilson(vx, vN, vCL=0.95):
	a = 1-vCL
	p = 1 - (a/2.)
	z = norm.ppf(p)
	phat = 1. * vx / vN
	X = z * sqrt((phat * (1 - phat) / vN) + (z * z / 4 / vN/ vN))
	invpre = 1 + (z * z / vN)
	brak_pre = phat + (z * z / (2 * vN))
	return (brak_pre - X) / invpre, (brak_pre + X) / invpre


if __name__ == '__main__':
	a,b,c = 162., 77., 4.
	import numpy as np
	N = np.sum([a,b,c])
	percentages = np.array([a,b,c]) / N

	confs = [wilson(i, N) for i in (a,b,c)]
	confs = [(percentages[ind] - i[0], i[1] - percentages[ind]) for ind, i in enumerate(confs)]
	
	for i, C in enumerate(confs):
		print " %s: %.2f + %.2f - %.2f" % (str(i), percentages[i], confs[i][1], confs[i][0])	