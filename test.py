import storage as S
import fit as F
from multiprocessing import Pool 

if __name__ == '__main__':
	tables, header = S.import_directory()
	N_list = range(3)

	images = [(tables[N], header[N]) for N in N_list]

	pool = Pool()
	pool.map(F.fit_bulge_disc, images)
	pool.close()
	pool.join()