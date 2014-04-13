import numpy as np
import time

map_time = []
for_time = []
# over_time = []
means = lambda x: np.mean(x)
for i in range(100):
	A = np.random.rand(1000, 50).tolist()
	s = time.clock()
	map(means, A)
	map_time.append(time.clock() - s)

	s = time.clock()
	[means(i) for i in A]
	for_time.append(time.clock() - s)

	# s = time.clock()
	# 1. / np.array(A)
	# over_time.append(time.clock() - s)

print 'map_time = %.2e +/- %.2e' % (np.mean(map_time), np.std(map_time))
print 'for_time = %.2e +/- %.2e' % (np.mean(for_time), np.std(for_time))
# print 'over_time = %.2e +/- %.2e' % (np.mean(over_time), np.std(over_time))

