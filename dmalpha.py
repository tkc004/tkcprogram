# it is the code to 

from readsnap_samson import *
from Sasha_functions import *
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from gadget_lib.cosmo import *
from samson_functions import *
mpl.use('Agg')
from matplotlib import rcParams
from dmalpha_tools import *
rcParams['axes.linewidth'] = 2
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['axes.unicode_minus']=False
rcParams.update({'figure.autolayout': True})

programdir='/home/tkc004/' #modify the directory to where tkcprogram locates

tlist=[13.9] #the cosmic time we consider
dirnr=[]
dirhr=[]
dirnr=['fm11v'] #the halo we consider here 
#dirnr = ['fm10q','fm11q','fm11v','f1146','f573','f553','f476','fm11','f383','f61','fm12c','fm12f','fm12i','fm12m','fm12q']
#dirhr = ['f1146_hv','f573_hv','f383_hv']
dirneed = np.concatenate((dirnr,dirhr),axis=0)
todo='alpha37' #alpha measured between 0.3kpc to 0.7kpc 
#todo='Malpha1251' #alpha measured between 0.01 to 0.02Rv
xaxis='ms_mv'
outputlist=1 #output a txt file for info
#plotname=todo+'_'+xaxis+'_'+'fire2'
plotname=todo+'_'+xaxis+'_'+'test' #name of the plot

if todo=='alpha37':
	alpha37list=[]
	mvirlist=[]
	mslist=[]
        alpha37list_hv=[]
        mvirlist_hv=[]
        mslist_hv=[]
	for ncount in range(len(dirneed)):
		runtodo = dirneed[ncount]
                haloinfo=cosmichalo(runtodo)
                highres=haloinfo['highres']
                if highres==1:
                        for icount in range(len(tlist)):
                                time=tlist[icount]
                                print 'time', time
                                Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.3,0.7],timesRv=0)
                                alpha37list_hv.append(Rvir37alpha)
                                mvirlist_hv.append(mvir)
                                mslist_hv.append(ms)
		else:
			for icount in range(len(tlist)):
				time=tlist[icount]
				print 'time', time
				Rvir37alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.3,0.7],timesRv=0)
				alpha37list.append(Rvir37alpha)
				mvirlist.append(mvir)
				mslist.append(ms)
	print 'dirneed', dirneed	
	print 'Mvir', mvirlist
	print 'Mstar', mslist
	print 'log10(Mstar)', np.log10(mslist)
	print 'alpha37list', alpha37list
	plotvsTHINGS(programdir)
	if len(alpha37list)>0:
		xneed=np.log10(mslist)
		yneed=alpha37list
		plotalpha37(xneed,yneed)
	if len(alpha37list_hv)>0:
		xneed=np.log10(mslist_hv)
		yneed=alpha37list_hv
		plotalpha37(xneed,yneed,marker='s',needlabel=0)
        finname = 'figures/'+plotname+'.pdf'
        plt.savefig(finname)
        plt.clf()


if todo=='Malpha1251':
	alpha51list=[]
        alpha12list=[]
        mvirlist=[]
        mslist=[]
        alpha51list_hv=[]
        alpha12list_hv=[]
        mvirlist_hv=[]
        mslist_hv=[]
        for ncount in range(len(dirneed)):
                runtodo = dirneed[ncount]
		haloinfo=cosmichalo(runtodo)
		highres=haloinfo['highres']
		if highres==1:
                        for icount in range(len(tlist)):
                                time=tlist[icount]
                                Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.01,0.02],timesRv=1)
                                Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.005,0.01],timesRv=1)
                                alpha12list_hv.append(Rvir12alpha)
                                alpha51list_hv.append(Rvir51alpha)
                                mvirlist_hv.append(mvir)
                                mslist_hv.append(ms)
		else:
			for icount in range(len(tlist)):
				time=tlist[icount]
				Rvir12alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.01,0.02],timesRv=1)
				Rvir51alpha, xdata, ydata, mvir, ms=dmalpha(runtodo,time,programdir,Rrange=[0.005,0.01],timesRv=1)
				alpha12list.append(Rvir12alpha)
				alpha51list.append(Rvir51alpha)
				mvirlist.append(mvir)
				mslist.append(ms)
	mslist=np.array(mslist)
	mvirlist=np.array(mvirlist)
        mslist_hv=np.array(mslist_hv)
        mvirlist_hv=np.array(mvirlist_hv)
        print 'dirneed', dirneed
        print 'Mvir', mvirlist
        print 'Mstar', mslist
        print 'log10(Mstar)', np.log10(mslist)
        print 'alpha12list', alpha12list
	print 'alpha51list', alpha51list
	if xaxis=='mv':
		if len(alpha12list)>0:
			xneed=mvirlist
			yneed=alpha12list
			y1need=alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis)
		if len(alpha12list_hv)>0:
			xneed=mvirlist_hv
			yneed=alpha12list_hv
			y1need=alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='s',needlabel=0)
	elif xaxis=='ms':
                if len(alpha12list)>0:
			xneed=mslist
			yneed=alpha12list
			y1need = alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis)
                if len(alpha12list)>0:
                        xneed=mslist_hv
                        yneed=alpha12list_hv
                        y1need = alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='s',needlabel=0)
	elif xaxis=='ms_mv':
                if len(alpha12list)>0:
			xneed=mslist/mvirlist
			yneed=alpha12list
			y1need=alpha51list
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis)
                if len(alpha12list)>0:
                        xneed=mslist_hv/mvirlist_hv
                        yneed=alpha12list_hv
                        y1need=alpha51list_hv
			plotMvMsRvalpha(xneed,yneed,y1need,zr=0.,xaxis=xaxis,marker='s',needlabel=0)
        finname = 'figures/'+plotname+'.pdf'
        plt.savefig(finname)
        plt.clf()


if outputlist==1:
        if todo=='Malpha1251':
		f = open('outputdata/dmalpha12Rv_info.txt', 'w')
		f.write('halo name    '+ 'Mvir    '+   'Mstar    '+ 'Mstar/Mvir   '+ 'alpha (0.01-0.02Rvir) '+ 'alpha (0.005-0.01) '+'\n')
		for i in range(len(mslist)):
			f.write(str(dirnr[i])+'      '+str(mvirlist[i])+'      '+str(mslist[i])+'      '+str(mslist[i]/mvirlist[i])+'      '+str(alpha12list[i])+'      '+str(alpha51list[i])+'  \n')
		for i in range(len(mslist_hv)):
			f.write(str(dirhr[i])+'      '+str(mvirlist_hv[i])+'      '+str(mslist_hv[i])+'      '+str(mslist_hv[i]/mvirlist_hv[i])+'      '+str(alpha12list_hv[i])+'      '+str(alpha51list_hv[i])+ '  \n')
        if todo=='alpha37':
		f = open('outputdata/dmalpha37_info.txt', 'w')
		f.write('halo name    '+ 'Mvir    '+   'Mstar    '+ 'Mstar/Mvir   '+ 'alpha (0.3-0.7kpc)'+'\n')
                for i in range(len(mslist)):
                        f.write(str(dirnr[i])+'      '+str(mvirlist[i])+'      '+str(mslist[i])+'      '+str(mslist[i]/mvirlist[i])+'      '+str(alpha37list[i])+'      '+'  \n')
                for i in range(len(mslist_hv)):
                        f.write(str(dirhr[i])+'      '+str(mvirlist_hv[i])+'      '+str(mslist_hv[i])+'      '+str(mslist_hv[i]/mvirlist_hv[i])+'      '+str(alpha37list_hv[i])+'      '+'  \n')
