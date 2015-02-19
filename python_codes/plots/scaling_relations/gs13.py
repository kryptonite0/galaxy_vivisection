import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy

def mbh_vs_kmag_sph():
	
	dataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/GS13/gs13.dat' 
	data = open(dataFileName)
	
	gal_nameList = []
	coreList = []
	mbhList = []
	perr_mbhList = []
	merr_mbhList = []
	kmag_sphList = []
	TypeList = []
	
	for line in data:
		if line.split()[0] != '#':
			gal_nameList.append(line.split()[0])
			coreList.append(line.split()[1])
			mbhList.append(10**8*float(line.split()[2]))
			perr_mbhList.append(10**8*float(line.split()[3]))
			merr_mbhList.append(10**8*float(line.split()[4]))
			kmag_sphList.append(float(line.split()[5]))
			TypeList.append(line.split()[6])
	
	data.close()
	
	gal_name = np.asarray(gal_nameList)
	core = np.asarray(coreList)
	mbh = np.asarray(mbhList)
	perr_mbh = np.asarray(perr_mbhList)
	merr_mbh = np.asarray(merr_mbhList)
	kmag_sph = np.asarray(kmag_sphList)
	Type  = np.asarray(TypeList)
	err_kmag_sph = kmag_sph*[0.0] 
	err_kmag_sph[Type=='E'] = +0.25
	err_kmag_sph[Type=='S'] = +0.75
	
	log_mbh = np.log10(mbh)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	
	fig, ax = plt.subplots()


	### fit using BCES ####
        print 'core-Sersic'
        A,B,Aerr,Berr,covAB=bces.bces(kmag_sph[core=='y']+25, err_kmag_sph[core=='y'], log_mbh[core=='y'],
        	0.5*(perr_log_mbh[core=='y'] + merr_log_mbh[core=='y']),kmag_sph[core=='y']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> = -25'
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
	
        logxx = np.arange(-28,-18,0.1)
	y_bcesbisec = (A[2]*(logxx+25) + B[2])
	ax.plot(logxx,y_bcesbisec, color='red', linewidth=2.)
	print logxx,y_bcesbisec
	
	### fit using BCES ####
        print 'Sersic'
        A,B,Aerr,Berr,covAB=bces.bces(kmag_sph[core=='n']+22.5, err_kmag_sph[core=='n'], log_mbh[core=='n'],
        	0.5*(perr_log_mbh[core=='n'] + merr_log_mbh[core=='n']),kmag_sph[core=='n']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> = -22.5'
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
	
	y_bcesbisec = (A[2]*(logxx+22.5) + B[2])
	ax.plot(logxx,y_bcesbisec, color='blue', linewidth=2.)
	
	### fit using FITEXY ###
        print 'core-Sersic'
        A,B = fitexy.bisect_modfitexy(kmag_sph[core=='y']+25, err_kmag_sph[core=='y'],
        	log_mbh[core=='y'], 0.5*(perr_log_mbh[core=='y'] + merr_log_mbh[core=='y']))
        # plot bisectore relation
        y_bisec = A[2] + B[2]*(logxx+25)
        ax.plot(logxx, y_bisec, ls='--', color='red', linewidth=2.)
	slope_cS = B[2]
      
        print 'Sersic'
        A,B = fitexy.bisect_modfitexy(kmag_sph[core=='n']+22.5, err_kmag_sph[core=='n'],
        	log_mbh[core=='n'], 0.5*(perr_log_mbh[core=='n'] + merr_log_mbh[core=='n']))
        # plot bisectore relation
        y_bisec = A[2] + B[2]*(logxx+22.5)
        ax.plot(logxx, y_bisec, ls='--', color='blue', linewidth=2.)
	slope_S = B[2]
	print 'slope_cS=', slope_cS
	print 'slope_S=', slope_S

	ax.errorbar(kmag_sph[core=='y'], log_mbh[core=='y'], xerr=err_kmag_sph[core=='y'], yerr=[perr_log_mbh[core=='y'],merr_log_mbh[core=='y']], ecolor='gray', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
	ax.errorbar(kmag_sph[core=='n'], log_mbh[core=='n'], xerr=err_kmag_sph[core=='n'], yerr=[perr_log_mbh[core=='n'],merr_log_mbh[core=='n']], ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
	ax.axis([-18,-28,5,11])
	plt.xlabel(r'$KMAG$', labelpad=20)
	plt.ylabel(r'$log(Mbh)$', labelpad=20)
	
	#plt.show()
        plt.savefig('mbh_vs_kmag_sph_GS13.pdf', format='pdf', dpi=1000)
	
	
	
	
mbh_vs_kmag_sph()	
	
