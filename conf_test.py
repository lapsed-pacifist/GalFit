from scipy.stats import norm
from numpy import sqrt
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman'], 'size':17})
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
	f = vx / vN
	a = 1-vCL
	p = 1 - (a/2.)
	z = norm.ppf(p)
	phat = 1. * vx / vN
	X = z * sqrt((phat * (1 - phat) / vN) + (z * z / 4 / vN/ vN))
	invpre = 1 + (z * z / vN)
	brak_pre = phat + (z * z / (2 * vN))
	return abs(f - ((brak_pre - X) / invpre)), abs(((brak_pre + X) / invpre) - f)


if __name__ == '__main__':
	import numpy as np
	erwin = [0.458*24., 0., (1-0.458)*24.]
	erwinN = 24.
	fieldN = 55.
	mineN = 66.
	mine = [16., 1., 49.]
	field = [0.26*55., 0.28*55., (1.-(0.28+0.26))*55.]
	mine_errs = [wilson(i, mineN, 0.68) for i in mine]
	erwin_errs = [wilson(i, erwinN, 0.68) for i in erwin]
	field_errs = [wilson(i, fieldN, 0.68) for i in field]
	
	fig = plt.figure()
	fig.set_facecolor('white')
	ax = fig.add_subplot(111)
	ind = np.arange(3)
	width = 0.8
	rects1 = ax.bar(ind, np.array(mine)/mineN, width/3., color='0.8', yerr=zip(*mine_errs), ecolor='0.1', label='Coma (This Study)')
	rects2 = ax.bar(ind + (width/3.), np.array(erwin)/erwinN, width/3., color='0.8', yerr=zip(*erwin_errs), ecolor='0.1', label='Erwin Virgo')
	rects3 = ax.bar(ind+(2*width/3.), np.array(field)/fieldN, width/3., color='0.4', yerr=zip(*field_errs), ecolor='0.1', label='Erwin Field')
	for b in rects2:
		b.set_hatch('\\')
	ax.set_xticks(ind+(width/2.))
	ax.set_xticklabels( ('Type I', 'Type II', 'Type III') )
	offset = 0.2
	ax.set_xlim([ax.get_xlim()[0]-offset, ax.get_xlim()[1]])
	ax.set_ylim([0,0.8])
	ax.set_ylabel('Fraction of S0')
	ax.legend(loc='best', prop={'size':16})
	print mine_errs
	plt.show()
