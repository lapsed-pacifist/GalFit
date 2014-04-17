import sys
import time
from numpy import floor

class loadbar():
	"""initiate with .time() at the start of each loop, do .progress at the end"""
	def __init__(self, length):
		self.length = length
		self.ind = 0
		if self.length < 2:
			self.nan = True
			print "Loadbar not viable for length <2"
		else:
			self.nan = False

	def time(self):
		if self.nan: return
		self.start_time = time.clock()

	def time_type(self, t):
		h = t/60/60
		m = (h - floor(h)) * 60
		s = (m - floor(m)) * 60
		h, m, s = map(floor, [h,m,s])
		if h == 0:
			if m == 0:
				return "%.f secs" % (s)
			else:
				return "%.f mins %.f secs" % (m, s)
		else:
			return "%.f:%.f:%.f hrs" % (h, m, s)


	def progress(self):
		if self.nan: return
		if self.ind == 0:
			self.av_t = time.clock() - self.start_time
			sys.stdout.write('\r')
			sys.stdout.write('First Iteration: Calculating...')
		else:
			 self.av_t = ((self.av_t*self.ind) + time.clock() - self.start_time) / (self.ind+1.)
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s] %.1f%%" % ('='*int(20*self.ind/(self.length-1)), (100*self.ind/(self.length-1))))
		sys.stdout.write(" left: %s" % (self.time_type(self.av_t*(self.length-self.ind))))
		sys.stdout.write(" at %.1f sec/iter (i=%i)    " % (self.av_t, self.ind+1))
		sys.stdout.flush()
		self.ind += 1

if __name__ == '__main__':
	select = range(10)

	S = time.clock()
	L = loadbar(len(select))
	for ind, i in enumerate(select):
		L.time()
		# >>do stuff
		time.sleep(1)
		L.progress()




