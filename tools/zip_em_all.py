import numpy as np

def zip_em_all(biglist, nfiles):
	if (nfiles < 2):
		print 'returning 0th element'
		zippedlist = np.array(biglist[0])
		return zippedlist
	zippedlist = np.append(biglist[0], biglist[1], axis = 0)
	count = 0
	endgame = nfiles - 2
	while (count < endgame):
		zippedlist = np.append(zippedlist, biglist[count+2], axis =0)
		count += 1
	return zippedlist
