import numpy as np
from tasz import tfora
import sys

#programdir='/home/tkc004'
def convertTemp(Tb, Neb, rho, h):
	#converts temperature and density to physical units
	H_MASSFRAC = 0.76
	PROTONMASS  = 1.6726e-24
	GAMMA = 5.0/3.0
	BOLTZMANN  = 1.3806e-16
	
	UnitLength_in_cm=3.085678e21 / h
	UnitMass_in_g=1.989e43 / h
	UnitVelocity_in_cm_per_s=1.e5
	UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
	UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)

	MeanWeights = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) * PROTONMASS
	MeanWeights_con = 4.0/(3*H_MASSFRAC+1+4*H_MASSFRAC* Neb) 

	converted_rho = rho * UnitMass_in_g * 6.02e23 / MeanWeights_con
	converted_rho /= (UnitLength_in_cm * UnitLength_in_cm *UnitLength_in_cm)
	
	TrueTemp = Tb *  MeanWeights/BOLTZMANN * (GAMMA-1.) * UnitEnergy_in_cgs/ UnitMass_in_g
	return (TrueTemp, converted_rho)

def cut_stellar_age(a, N, age, omeganot, littleh):
	#returns the time since the previous snapshot, as well as an array that gives the actual age of all stellar particles in Myr (as opposed to formation time, as given by code)
	finname = 'output_times.txt'
	f = open(finname)
	SdarsA = np.loadtxt(f)
	prev_epoch = SdarsA[max((N-1), 0)]
	SFR_time_range = (tfora(a, omeganot, littleh) - tfora(prev_epoch, omeganot, littleh)) * 1e9
	print 'using exact SFR', a, prev_epoch, (SFR_time_range)/1e7, ' * 1e7 yrs'
	f.close()
	#print age
	time_in_Gyr = tfora(age, omeganot, littleh)
	age_in_Gyr = tfora(a, omeganot, littleh) - time_in_Gyr
	print ' i interpolated your ages'
	#print age_in_Gyr
	#print ' the times of formation are: ',time_in_Gyr
	return (SFR_time_range, age_in_Gyr)


def read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',maindir='scratch'):
	finname = maindir+'/'+rundir+'/Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_halos' #or whatever prefix
	f = open(finname)
	C = np.loadtxt(f)
	halID = C[:,0]
	x = C[:,5]
	y = C[:,6]
	z = C[:,7]
	vx = C[:,8]
	vy = C[:,9]
	vz = C[:,10]
	M = C[:,3]
	Vmax = C[:,16]
	Vsig = C[:,18]
	Rvir = C[:,11]
	Mgas = C[:,53]
	Mstar = C[:,73]
	fMhires = C[:,37]
	#if (use_mbp=='y'): not even going to bother writing this for now - it isnt very accurate for halo center in Amiga. May want to use COM instead though
		
	return {'fMhires':fMhires,'ID':halID, 'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'M':M, 'Vmax':Vmax, 'Vsig':Vsig, 'Rvir':Rvir, 'Mgas':Mgas, 'Mstar':Mstar}

	
def shell_check(dists, Rvir, R_bins, a, h, use_physical_bins='n', physical_extent=50):

	TheCuts = []
	bincount = 0
	ShellList = np.zeros(len(dists)) - 1 
	while (bincount < R_bins):
		DistMin = Rvir * bincount/R_bins
		DistMax = Rvir * (bincount+1)/R_bins
	
		#by default, the virial radius is cut into 10 bins. Here i'm also allowing for a fixed physical scale.
		if (use_physical_bins == 'y'):
			print 'using physical bins ',physical_extent, 'kpc/h'
			#gotta divide by a to go to comoving coordinates
			Distmin = bincount/(Rbins) * physical_extent / a
			Distmax = (bincount+1)/(Rbins) * physical_extent / a 

		ShellThickness = (DistMax - DistMin)*a/h

		Cut1 = dists > DistMin
		Cut2 = dists < DistMax
		Cuts = Cut1 * Cut2
		ShellList[Cuts] = bincount
		TheCuts.append(Cuts)
		bincount+=1
	TheCuts = np.array(TheCuts)
	return (TheCuts, ShellList, ShellThickness)

def calcdist2(x1, y1, z1, posx, posy, posz, boxsize):
#calculates distances of an array of points (particles) from a single given point (halo center)
#for your reference: x1, y1, z1 are floats for halo center or hwatever, posx, posy posz are arrays for positions of particles or whatever
	 xdist = abs(posx - x1)
	 ydist = abs(posy - y1)
	 zdist = abs(posz - z1)
	 
	 #adjust for periodic boundary conditions
	 x_straggler = xdist/boxsize > 0.5
	 if (len(xdist[x_straggler])) > 0:
	 	print 'adjusting x stragglers ', len(xdist[x_straggler])
	 	xdist[x_straggler] = boxsize - xdist[x_straggler]
	 y_straggler = ydist/boxsize > 0.5
	 if (len(ydist[y_straggler])) > 0:
	 	ydist[y_straggler] = boxsize - ydist[y_straggler]	 	
	 	print 'adjusting y stragglers ', len(ydist[y_straggler])
	 z_straggler = zdist/boxsize > 0.5
	 if (len(zdist[z_straggler])) > 0:
	 	zdist[z_straggler] = boxsize - zdist[z_straggler]
	 	print 'adjusting z stragglers ', len(zdist[z_straggler])

	 
	 dists = pow(((xdist*xdist + ydist*ydist + zdist*zdist)), 0.5)
 	 return dists
 	 
def line_to_string(line, newline = 'y'):
	linestring = ''
	for q in line:
		linestring += "{0:.6g}".format(q) + '  '
	if (newline == 'y'):
		linestring+=' \n'
	return linestring
	
def sasha_max(a):
	if (len(a) > 0):
		return np.max(a)
	else:
		print 'warning there was an empty array max'
		return -99999999

def sasha_min(a):
	if (len(a) > 0):
		return np.min(a)
	else:
		print 'warning there was an empty array min'
		return 99999999

def check_Rvir_growth(halocount, a, Rvir, Vsig, Mass):
	#check to make sure that the snapshot you are using is not affected by a glitch in the halo catalog that causes a temporary dip in Rvir, Mvir, Vsig 
	H = read_halo_history(int(halocount))
	history_Rvirs = H['Rvir']
	history_M = H['M']
	history_Vsig = H['Vsig']
	history_Zs = H['redshift']
	history_As = 1.0 / (1.0 + history_Zs)
	history_Rvirs_phys = history_Rvirs * history_As
	RvirPhys = Rvir * a
	cut = history_As < a
	historical_max = sasha_max(history_Rvirs_phys[cut])  
	if (RvirPhys >= historical_max):
		print 'Rvir is already the biggest'
		return (Rvir, Vsig, Mass)
	maxindex = np.where(history_Rvirs_phys[cut]==historical_max)[0][0]
	print 'adjusting back to z = ',history_Zs[maxindex]
	print 'Rvir was ',history_Rvirs[maxindex]
	print 'instead of meager ',Rvir
	return (history_Rvirs[maxindex], history_Vsig[maxindex], history_M[maxindex])

def read_halo_history(rundir, halonostr='00', multifile='n', hubble=1,comoving=1, maindir='scratch'): #to output comoving distance comoving = 1
        redlist = []
        halolist = []
        xlist = []
        ylist = []
        zlist = []
        mvirlist = []
	mstarlist = []
	mgaslist = []
	rvirlist = []
        fMhireslist = []
	dirname='/home/tkc004/'+maindir+'/'+rundir
	halofile=open(dirname+'/halos/halo_000'+halonostr+'.dat','r')
	halofile.readline()
	halofile.readline()
	dars = halofile.readlines()
	halofile.close()
	for line in dars:
		xsd = line.split()
		mvir = float(xsd[4])
		xcen = float(xsd[6])
		ycen = float(xsd[7])
		zcen = float(xsd[8])
		if comoving==1:
			zred=0.0
		else:
			zred=float(xsd[0])
		redlist=np.append(redlist,float(xsd[0]))
		halolist=np.append(halolist,int(xsd[1]))
		mvirlist=np.append(mvirlist,mvir/hubble)
		mstarlist=np.append(mstarlist,float(xsd[74])/hubble)
		xlist=np.append(xlist,xcen/hubble/(1.0+zred)) # xcen originally in comoving unit (kpc/h)
		ylist=np.append(ylist,ycen/hubble/(1.0+zred))
		zlist=np.append(zlist,zcen/hubble/(1.0+zred))
		rvirlist=np.append(rvirlist,float(xsd[12])/hubble/(1.0+zred)) # originally in comoving unit (kpc/h)
                fMhireslist= np.append(fMhireslist, float(xsd[38]))
		mgaslist = np.append(mgaslist, float(xsd[54])/hubble)
	redlist=np.array(redlist)
	halolist=np.array(halolist)
	mvirlist=np.array(mvirlist)
	mstarlist=np.array(mstarlist)
	mgaslist=np.array(mgaslist)
	xlist=np.array(xlist)
	ylist=np.array(ylist)
	zlist=np.array(zlist)
	rvirlist=np.array(rvirlist)
        fMhireslist=np.array(fMhireslist)
	return {'redshift':redlist, 'ID':halolist, 'x':xlist, 'y':ylist, 'z':zlist, 'M':mvirlist, 'Ms':mstarlist, 'Mg':mgaslist, 'R':rvirlist, 'fMhires':fMhireslist}


def read_halo_history_pep(rundir, finalno, programdir,singlesnap=1, beginno=100,halonostr='00', multifile='n', hubble=1,comoving=1, maindir='scratch',firever=1): #to output comoving distance comoving = 1
        redlist = []
        idlist = []
        xlist = []
        ylist = []
        zlist = []
        Mlist = []
        Mslist = []
        Mglist = []
        Rvlist = []
        fMlist = []
	halono=int(halonostr)
	strnolist, zstrlist = snapshotzstr(programdir,firever=firever)
	if (singlesnap==1): beginno=finalno
	for Nsnap in range(beginno,finalno+1):
		redshiftstring = str(zstrlist[Nsnap])
		Nsnapstring = str(strnolist[Nsnap]) 
		hcat=read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',maindir=maindir)
		mvirl = hcat['M']
		xl = hcat['x']
		yl = hcat['y']
		zl = hcat['z']
		Msl = hcat['Mstar']
		Mgl = hcat['Mgas']
		Rl = hcat['Rvir']
		fMl = hcat['fMhires']
		redshift = float(redshiftstring)
		if comoving==1:
			zred=0.0
		else:
			zred=redshift
		haloid = halono
		mvir = float(mvirl[halono]/hubble)
		mstar = float(Msl[halono]/hubble)
		mgas = float(Mgl[halono]/hubble)
		xcen = float(xl[halono]/hubble/(1.0+zred)) # xcen originally in comoving unit (kpc/h)
		ycen = float(yl[halono]/hubble/(1.0+zred))
		zcen = float(zl[halono]/hubble/(1.0+zred))
		rvir = float(Rl[halono]/hubble/(1.0+zred))  # originally in comoving unit (kpc/h)
		xlist = np.append(xlist,xcen)
		ylist = np.append(ylist,ycen)
		zlist = np.append(zlist,zcen)
		idlist = np.append(idlist,haloid)
		redlist = np.append(redlist,float(zstrlist[Nsnap]))
		Mlist = np.append(Mlist,mvir)
		Mslist = np.append(Mslist,mstar)
		Mglist = np.append(Mglist,mgas)
		Rvlist = np.append(Rvlist,rvir)
		fMlist = np.append(fMlist,fMl)
	xlist = np.array(xlist)
	ylist = np.array(ylist)
	zlist = np.array(zlist)
	idlist = np.array(idlist)
	redlist = np.array(redlist)
	Mlist = np.array(Mlist)
	Mslist = np.array(Mslist)
	Mglist = np.array(Mglist)
	Rvlist = np.array(Rvlist)
	fMlist = np.array(fMlist)
	if (singlesnap==1):
		return {'fMhires':fMl, 'redshift':redshift, 'ID':haloid, 'x':xcen, 'y':ycen, 'z':zcen, 'M':mvir, 'Ms':mstar, 'Mg':mgas, 'R':rvir}
	else:
		return {'fMhires':fMlist, 'redshift':redlist, 'ID':idlist, 'x':xlist, 'y':ylist, 'z':zlist, 'M':Mlist, 'Ms':Mslist, 'Mg':Mglist, 'R':Rvlist}
def snapshotzstr(programdir,firever=1):
	if firever==1:
		print 'FIRE 1 under construction'
		exit()
	spmname=programdir+'/tkcprogram/data/snapshot_times.txt'
	spmfile=open(spmname,"r")
	spmfile.readline()
	spmfile.readline()
	spmfile.readline()
	dars = spmfile.readlines()

	snapno=[]
	redshift=[]

	for line in dars:
		xsd = line.split()
		snapno = np.append(snapno, int(xsd[0]))
		redshift = np.append(redshift, float(xsd[2]))
	spmfile.close()
	snapno=np.array(snapno)
	redshift=np.array(redshift)
	strzlist=[]
	strnolist=[]
	for icount in range(len(snapno)):
		if snapno[icount]>-1:
			strno=str(int(snapno[icount]))
			if icount <10:
				strno =  '00' + strno
			elif icount<100:
				strno = '0'+strno
			if strno=='250':
				strz = '1.187'
			elif strno=='000':
				strz = '99.013'
			elif strno=='024':
				strz = '9.313'
			else:
				strz="%0.3f" % redshift[icount]
		strzlist = np.append(strzlist,strz)
		strnolist = np.append(strnolist,strno)
	return strnolist, strzlist
