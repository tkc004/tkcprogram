from struct import *
import sys
import os
import matplotlib.pyplot as plt
import pylab
import numpy as np
import math
from distcalcs import *
from zip_em_all import *
import time
from scipy.interpolate import interp1d
from datetime import date
import scipy.stats as stats
import scipy.optimize as optimize
import errno
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
from Sasha_functions import read_halo_history
from readsnap_samson import *

#programdir='/home/tkc004/' #modify the directory to where tkcprogram locates

pi  = np.pi
sin = np.sin
cos = np.cos

def ellipse(u,v):
    x = rx*cos(u)*cos(v)
    y = ry*sin(u)*cos(v)
    z = rz*sin(v)
    return x,y,z


def mvee(points, tol = 0.001):
    """
    Finds the ellipse equation in "center form"
    (x-c).T * A * (x-c) = 1
    """
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T
    err = tol+1.0
    u = np.ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = np.dot(np.dot(Q, np.diag(u)), Q.T)
        M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
    c = np.dot(u,points)
    A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points)
               - np.multiply.outer(c,c))/d
    return A, c


# This routine uses relative positions of the particles to locate the center. The basic mechanism is to group particles into rectangular grids, and then fit the high density grid with ellipsoidal. 

#DMrelpos is the positions of the particles and nobin^3 is the number of grid. 

#nopafactor* den is the density cutoff (den is the average density of the particles). 

# !It is better to keep nobin <100 or the code will cost a lot of time and the center will be locating the highest density clump but not based on overall shape 

# ! It is not possible to fit with too many or too few grids (if nopa is too small, the number of grid will be too large and you should tune the nopafactor


def ellipsoidal_centering(DMrelposX,DMrelposY,DMrelposZ,nobin,nopafactor):
	DMco=np.array([DMrelposX,DMrelposY,DMrelposZ])
	DMcor=DMco.T

	DMparno=len(DMrelposX)
	H, edges = np.histogramdd(DMcor, bins = (nobin, nobin, nobin))

	NominS=[]
	ellX=[]
	ellY=[]
	ellZ=[]

	den=float(DMparno)/(pi*4./3.)/float(nobin)/float(nobin)/float(nobin)

	nopa=float(den)*nopafactor

	Dpos=[]
	Dposx=[]
	Dposy=[]
	Dposz=[]
	totalno=0
	for i in range(nobin):
		for j in range(nobin):
			for k in range(nobin):
				if (H[i,j,k]>nopa):
					Dposx=np.append(Dposx,edges[0][i])
					Dposy=np.append(Dposy,edges[1][j])
					Dposz=np.append(Dposz,edges[2][k])
					totalno+=H[i,j,k]
	if len(Dposx)<4:
		print 'Warning: Density threshold is too high'
		return -1, -1, -1 
	if len(Dposx)>1000:
		print 'Warning: Density threshold is too low and the no of grids is too large'
		return -1, -1, -1

	print 'fitting ell'
	points = np.vstack(((Dposx),(Dposy),(Dposz))).T
	A, centroid = mvee(points)
	DXi=centroid[0]
	DYi=centroid[1]
	DZi=centroid[2]
	return DXi, DYi, DZi



# it is a code to read halo_000xx.dat generated with AHF and merger tree

def readhalos(dirname,halostr,hubble=0.702):
        halofile=open(dirname+'/halos/halo_000'+halostr+'.dat','r')
        halofile.readline()
        halofile.readline()
        dars = halofile.readlines()
        halofile.close()
        zlist=[]
        idlist=[]
        mvirlist=[]
        fMhireslist=[]
        mgaslist=[]
        mstarlist=[]
        rvirlist=[]
        haloXl=[]
        haloYl=[]
        haloZl=[]
        for line in dars:
                xsd = line.split()
                zlist.append(float(xsd[0]))
                idlist.append(int(xsd[1]))
                mvirlist.append(float(xsd[4])/hubble)
                haloXl.append(float(xsd[6])/hubble)
                haloYl.append(float(xsd[7])/hubble)
                haloZl.append(float(xsd[8])/hubble)
                fMhireslist.append(float(xsd[38]))
                mgaslist.append(float(xsd[54])/hubble)
                mstarlist.append(float(xsd[74])/hubble)
                rvirlist.append(float(xsd[12])/hubble)
	zlist=np.array(zlist)
	a_scale = 1./(1.+zlist)
	idlist=np.array(idlist)
	mvirlist=np.array(mvirlist)
	fMhireslist=np.array(fMhireslist)
	#print 'a_scale', a_scale
	#print 'haloXl', haloXl
	haloX=np.array(haloXl)*a_scale
	haloY=np.array(haloYl)*a_scale
	haloZ=np.array(haloZl)*a_scale
	#print 'haloX', haloX
	mgaslist=np.array(mgaslist)
	mstarlist=np.array(mstarlist)
	rvirlist=np.array(rvirlist)
	return {'k':1,'mv':mvirlist,'fM':fMhireslist,'haloX':haloX,'haloY':haloY,'haloZ':haloZ,'mg':mgaslist,'ms':mstarlist,'rv':rvirlist,'z':zlist,'id':idlist};	




#it is a code to locate the center of a halo/galaxy (SXi,SYi,SZi) and get the effective radius (rx,ry,rz) through an ellipsoid fit. 
#But it is different from ellipsoidal_centering, because it only counts grids above a certain density threshold such that the ellipse encloses a fraction of the total mass (upmr and dmr are the upper and lower bounds of the mass fraction), but ellipsoidal_centering gives the ellipse that encloses ALL of the particles.

#H is a 3 dimensional array that contains density information and edges are the spatial boundary of H. One way to generate them is to use np.histogramdd, e.g. see the function ellipsoidal_centering above

#haloX, haloY and haloZ are the initial guess of the center
#nobin is the size of H[i,j,k]
#den (times 200) is the initial guess of the density threshold

def reff_ell(H, edges, haloX, haloY, haloZ, nobin, den):
        NominS=[]
        ellX=[]
        ellY=[]
        ellZ=[]
        Spos=[]
        Sposx=[]
        Sposy=[]
        Sposz=[]
	errlist=[]
        SXi=haloX
        SYi=haloY
        SZi=haloZ
        nit=0
        massrat=0
        upno=1000
        upmr=0.45
        dmr=0.55
        nopa=float(den*200)
        NominS=np.append(NominS,nopa)
        nopau=float(den*400)
        nopad=float(den*10)
        inell=0
        while (nit<1000 and (massrat<dmr or massrat>upmr)):
                nit+=1
                if inell==1:
                        nopa=(nopau+nopad)/2.
                NominS=np.append(NominS,nopa)

                print 'nopa', nopa
                print 'len(Sposx)', len(Sposx)
                Spos=[]
                Sposx=[]
                Sposy=[]
                Sposz=[]
                totalno=0
                for i in range(nobin):
                        for j in range(nobin):
                                for k in range(nobin):
                                        if (H[i,j,k]>nopa):
                                                Sposx=np.append(Sposx,edges[0][i])
                                                Sposy=np.append(Sposy,edges[1][j])
                                                Sposz=np.append(Sposz,edges[2][k])
                                                totalno+=H[i,j,k]
                print 'len(Sposx)', len(Sposx)
                if ((len(Sposx)>3 and len(Sposx)<upno) or inell==1):
                        inell=1
                        #print 'fitting ell'
                        points = np.vstack(((Sposx+haloX),(Sposy+haloY),(Sposz+haloZ))).T
                        A, centroid = mvee(points)
                        print 'centroid', centroid
                        SXi=centroid[0]
                        SYi=centroid[1]
                        SZi=centroid[2]
                        U, D, V = la.svd(A)
                        rx, ry, rz = 1./np.sqrt(D)
                        u, v = np.mgrid[0:2*pi:20j, -pi/2:pi/2:10j]

                        def ellipse(u,v):
                            x = rx*cos(u)*cos(v)
                            y = ry*sin(u)*cos(v)
                            z = rz*sin(v)
                            return x,y,z

                        edgespoints=[]
                        E = np.dstack(ellipse(u,v))
                        E = np.dot(E,V) + centroid
                        x, y, z = np.rollaxis(E, axis = -1)
                        err=0
                        errlist=np.append(errlist,0)
                        inV=la.inv(V)
                        for i in range(1,len(edges[0])):
                                for j in range(1,len(edges[1])):
                                        for k in range(1, len(edges[2])):
                                                edgespoints=np.append(edgespoints,(((edges[0][i]+edges[0][i-1])/2-centroid[0]+haloX),((edges[1][j]+edges[1][j-1])/2-centroid[1]+haloY),((edges[2][k]+edges[2][k-1])/2-centroid[2]+haloZ)))
                        edgespoints=np.matrix(edgespoints.reshape((len(edgespoints)/3,3)))
                        rotback=np.dot(edgespoints,inV).T
                        sumHcrit=0
                        for i in range(0,len(edges[0])-1):
                                for j in range(0,len(edges[1])-1):
                                        for k in range(0, len(edges[2])-1):
                                                n=i*(len(edges[2])-1)*(len(edges[2])-1)+j*(len(edges[2])-1)+k
                                                if (np.power(rotback[0,n]/rx,2)+np.power(rotback[1,n]/ry,2)+np.power(rotback[2,n]/rz,2)<1):
                                                        sumHcrit=sumHcrit+H[i,j,k]
                        sumH=np.sum(H)
                        massrat=sumHcrit/sumH
                        print 'sum in ell/all', massrat
                        del x, y, z, E, points, A, centroid,U, D, V,u, v, inV, rotback
                        massratm=massrat
                        nopam=nopa
                        if massratm>upmr:
                                nopad=nopa*1.2
                        if massratm<dmr:
                                nopau=nopa*0.8
                        if (nit>3):
#                               print 'NominS', NominS[nit-1], NominS[nit-3]
                                if (np.abs(NominS[nit-1]-NominS[nit-3])<0.0001*NominS[nit-1]):
                                        err=5
                                        errlist=np.append(errlist,5)
                                        break
                        if (sumHcrit<1e-8):
                                err=6
                                errlist=np.append(errlist,6)
                                break
                elif (inell==0):
                        if(len(Sposx)<3):
                                nopa=nopa*0.9-1
                                nopau=nopa
                        if(len(Sposx)>upno):
                                nopa=nopa*1.1+1
                                nopad=nopa
	return SXi, SYi, SZi, rx, ry, rz, err, massrat, Sposx


#it is a function that stores information of each simulation.


def cosmichalo(runtodo):
	multifile='n' #whether the hdf5 for one snapshot is a single or multiple files
        fileno=4 # if multifile='y', then we need the number of hdf5 files for one snapshot
	subdir='hdf5/' #the directory of hdf5 file (within the run directory)
	# complete path for hdf5 would be: maindir/rundir/subdir
	beginno=400 #the first snapshot number to consider
	finalno=441 # the last snapshot number to consider
	halocolor='k' #the color of the halo in plot
	labelname='not set' # the label of the halo in plot
	xmax_of_box=40.0  #only useful for UDG image
	halostr='00' # the halo number from AHF finder and Mergertree (e.g. halo_00000.dat -> '00')
	firever=1 # FIRE-1 -> 1; FIRE-2 -> 2
	usepep=0 #use AHF .halos instead of halo_0000x.dat
	maindir='/home/tkc004/scratch/' 
	highres=0 

        if (runtodo=='fm11v'):
                rundir='m11v/'
                halostr='00'
                beginno=100
                finalno=600
                subdir='output'
                xmax_of_box=40.0
                halocolor='b'
                labelname='fm11v'
                firever=2
                maindir='/home/tkc004/oasis/extra/'
                usepep=1


	return {'fileno':fileno,'rundir':rundir,'subdir':subdir,'halostr':halostr,'beginno':beginno,'finalno':finalno, 'multifile':multifile, 'halocolor':halocolor, 'labelname':labelname, 'xmax_of_box':xmax_of_box, 'firever':firever, 'usepep':usepep, 'maindir':maindir, 'highres':highres}


#to find the effective radius 

def reff_ell2d(H, edges, haloX, haloY, nobin, den):
        NominS=[]
        ellX=[]
        ellY=[]
        Spos=[]
        Sposx=[]
        Sposy=[]
        errlist=[]
        SXi=haloX
        SYi=haloY
        nit=0
        massrat=0
        upno=1000
        upmr=0.3
        dmr=0.7
        nopa=float(den)
        NominS=np.append(NominS,nopa)
        nopau=float(den*2)
        nopad=float(den*0.2)
        inell=0
	rx=0.
	ry=0.
	SXi=0.
	SYi=0.
        while (nit<100 and (massrat<dmr or massrat>upmr)):
                nit+=1
                if inell==1:
                        nopa=(nopau+nopad)/2.
                NominS=np.append(NominS,nopa)

                print 'nopa', nopa
                print 'len(Sposx)', len(Sposx)
                Spos=[]
                Sposx=[]
                Sposy=[]
                totalno=0
                for i in range(nobin):
                        for j in range(nobin):
				if (H[i,j]>nopa):
					Sposx=np.append(Sposx,edges[0,i])
					Sposy=np.append(Sposy,edges[1,j])
					totalno+=H[i,j]
                print 'len(Sposx)', len(Sposx)
 		print 'totalno', totalno
                if ((len(Sposx)>3 and len(Sposx)<upno) or inell==1):
                        inell=1
		 	#print 'Sposx', Sposx
                        #print 'fitting ell'
                        points = np.vstack(((Sposx+haloX),(Sposy+haloY))).T
                        A, centroid = mvee(points)
                        print 'centroid', centroid
                        SXi=centroid[0]
                        SYi=centroid[1]
                        U, D, V = la.svd(A)
                        rx, ry = 1./np.sqrt(D)
                        u = np.mgrid[0:2*pi:20j]

                        def ellipse2d(u):
                            x = rx*cos(u)
                            y = ry*sin(u)
                            return x,y

                        edgespoints=[]
                        E = np.dstack(ellipse2d(u))
                        E = np.dot(E,V) + centroid
                        x, y= np.rollaxis(E, axis = -1)
                        err=0
                        errlist=np.append(errlist,0)
                        inV=la.inv(V)
                        for i in range(1,nobin):
                                for j in range(1,nobin):
					edgespoints=np.append(edgespoints,(((edges[0,i]+edges[0,i-1])/2-centroid[0]+haloX),((edges[1,j]+edges[1,j-1])/2-centroid[1]+haloY)))
                        edgespoints=np.matrix(edgespoints.reshape((len(edgespoints)/2,2)))
                        rotback=np.dot(edgespoints,inV).T
                        sumHcrit=0
                        for i in range(0,len(edges[0])-1):
                                for j in range(0,len(edges[1])-1):
					n=i*(len(edges[1])-1)+j
					if (np.power(rotback[0,n]/rx,2)+np.power(rotback[1,n]/ry,2)<1):
						sumHcrit=sumHcrit+H[i,j]
                        sumH=np.sum(H)
                        massrat=sumHcrit/sumH
                 #       print 'sum in ell/all', massrat
                        del x, y, E, points, A, centroid,U, D, V,u, inV, rotback
                        massratm=massrat
                        nopam=nopa
                        if massratm>upmr:
                                nopad=nopa*1.1
                        if massratm<dmr:
                                nopau=nopa*0.9
                        if (nit>3):
                #                print 'NominS', NominS[nit-1], NominS[nit-3]
                                if (np.abs(NominS[nit-1]-NominS[nit-3])<0.0001*NominS[nit-1]):
                                        err=5
                                        errlist=np.append(errlist,5)
                                        break
                        if (sumHcrit<1e-8):
                                err=6
                                errlist=np.append(errlist,6)
                                break
                elif (inell==0):
                        if(len(Sposx)<3):
                                nopa=nopa*0.97-1
                                nopau=nopa
                        if(len(Sposx)>upno):
                                nopa=nopa*1.03+1
                                nopad=nopa
        return SXi, SYi, rx, ry,  err, massrat, Sposx



def readtime(programdir,firever=2):
	file=open(programdir+'/tkcprogram/data/snapshot_times.txt','r')
	file.readline()
	file.readline()
	file.readline()
	dars = file.readlines()
	file.close()
	snap2list=[]
	a2list=[]
	time2list=[]
	red2list=[]
	for line in dars:
		xsd = line.split()
		snap2list.append(int(xsd[0]))
		a2list.append(float(xsd[1]))
		red2list.append(float(xsd[2]))
		time2list.append(float(xsd[3]))
	snap2list=np.array(snap2list)
	time2list=np.array(time2list)
	if firever==1:
		file=open(programdir+'/tkcprogram/data/output_times.txt','r')
		dars = file.readlines()
		file.close()
                snaplist=[]
                alist=[]
                timelist=[]
                redlist=[]
		ncount=0
                for line in dars:
                        xsd = line.split()
                        alist.append(float(xsd[0]))
			snaplist.append(ncount)
			ncount+=1
		alist=np.array(alist)
		snaplist=np.array(snaplist)
		timelist=np.array(np.interp(alist,a2list,time2list))
		return snaplist, timelist
	if firever==2:
		return snap2list, time2list

