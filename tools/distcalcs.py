from numpy import *

def calcdist(x1, y1, z1, poses):
	boxsize = 10000.0 
	dists = []
	for poslist in poses:

	 xdist = abs(x1 - poslist[0])
	 ydist = abs(y1 - poslist[1])
	 zdist = abs(z1 - poslist[2])
	 if ( (xdist / boxsize)  > 0.5):
			xdist = min( abs(x1  - (poslist[0] + boxsize)), abs( poslist[0] - (x1 + boxsize)))
	 if ( (zdist / boxsize)  > 0.5):
			zdist = min( abs(z1  - (poslist[2] + boxsize)), abs( poslist[2] - (z1 + boxsize)))		 
	 if ( (ydist / boxsize)  > 0.5):
			ydist = min( abs(y1  - (poslist[1] + boxsize)), abs( poslist[1] - (y1 + boxsize)))
	 thedist = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
	 dists.append(thedist)
	 print x1, y1, z1, poslist, thedist
 	return dists

def calcdist2(x1, y1, z1, posx, posy, posz):
	 xdist = abs(posx[0:] - x1)
	 ydist = abs(posy[0:] - y1)
	 zdist = abs(posz[0:] - z1)
	 dists = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
 	 return dists


def calcdistE(x1, y1, z1, x2, y2, z2):
	boxsize = 10000.0 / 0.7
	xdist = abs(x1 - x2)
	ydist = abs(y1 - y2)
	zdist = abs(z1 - z2)
	if ( (xdist / boxsize)  > 0.5):
			xdist = min( abs(x1  - (x2 + boxsize)), abs( x2 - (x1 + boxsize)))
	if ( (zdist / boxsize)  > 0.5):
			zdist = min( abs(z1  - (z2 + boxsize)), abs( z2 - (z1 + boxsize)))		 
	if ( (ydist / boxsize)  > 0.5):
			ydist = min( abs(y1  - (y2 + boxsize)), abs( y2 - (y1 + boxsize)))
	thedist = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
	return thedist	

def mag(posx, posy, posz):
	themag = posx*posx + posy*posy + posz*posz + 1e-24
	themag = pow(themag, 0.5)
	return themag


def middlepointarray(oarray):
	last=len(oarray)-1
	sarray=delete(oarray, 0, axis=0)
	larray=delete(oarray, last, axis=0)
	tmp=(sarray+larray)/2
	narray=array(tmp)
	return narray

def diffarray(oarray):
        last=len(oarray)-1
        sarray=delete(oarray, 0, axis=0)
        larray=delete(oarray, last, axis=0)
        tmp=absolute(sarray-larray)
        narray=array(tmp)
        return narray

