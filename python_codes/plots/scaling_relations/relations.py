import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments import b_n

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 22})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_auxiliary_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/auxiliary/'
path_scalrel_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/scaling_relations/'


def mag_sph_vs_logn(axis):
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.ELLIPTICAL_my, \
		pysres.mag_sph_eq_moffat_comb, pysres.log_n_' + axis + '_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, errV.perr_log_n, errV.merr_log_n \
        	FROM Ancillary AS anc \
        	JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
        	WHERE anc.fit1D_done = 1;'
       
        cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
        #print data
        ELL = data[1].astype(np.int)
        mag_sph = data[2].astype(np.float)
        log_n = data[3].astype(np.float)
	perr_mag_sph = data[4].astype(np.float)
	merr_mag_sph = data[5].astype(np.float)
	perr_log_n = data[6].astype(np.float)
	merr_log_n = data[7].astype(np.float)
       
        fig, ax = plt.subplots()
        ax.errorbar(log_n[ELL==1], mag_sph[ELL==1], xerr=[merr_log_n[ELL==1],perr_log_n[ELL==1]], yerr=[merr_mag_sph[ELL==1],perr_mag_sph[ELL==1]], 
		ecolor='gray', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[ELL==0], mag_sph[ELL==0], xerr=[merr_log_n[ELL==0],perr_log_n[ELL==0]], yerr=[merr_mag_sph[ELL==0],perr_mag_sph[ELL==0]], 
		ecolor='gray', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 

        xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
	
	print 'Ellipticals'
	#print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
	A,B,Aerr,Berr,covAB=bces.bces(log_n[ELL==1]-np.average(log_n[ELL==1]),0.5*(perr_log_n[ELL==1]+merr_log_n[ELL==1]),mag_sph[ELL==1],
		0.5*(merr_mag_sph[ELL==1] + perr_mag_sph[ELL==1]),log_n[ELL==1]*[0.0])
	print '---------------------------------'
	print 'y = A*(x-<x>) + B '
	print '<x> =', np.average(log_n[ELL==1])
	print
	print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
	print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
	print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
	print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
	print '---------------------------------'
	
	xx = np.arange(0.001,30,0.1)
	logxx = np.log10(xx)
	#print A[2], B[2]
	# calculates the prediction bands for the given input arrays
	lpb68,upb68,logxx = predband.predband(log_n[ELL==1]-np.average(log_n[ELL==1]),mag_sph[ELL==1],A[2],B[2],conf=0.68,x=logxx)
	#lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
	#lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
	yy = (A[2]*(logxx) + B[2])
	
	ax.plot(logxx+np.average(log_n[ELL==1]),yy, color='red', linewidth=2.)
	# plots a shaded area containing the prediction band  
	plt.fill_between(logxx+np.average(log_n[ELL==1]), lpb68, upb68, alpha=0.15, facecolor='red')
	ax.plot(logxx+np.average(log_n[ELL==1]),lpb68, color='red')
	ax.plot(logxx+np.average(log_n[ELL==1]),upb68, color='red')

        print 'Bulges'
        A,B,Aerr,Berr,covAB=bces.bces(log_n[ELL==0]-np.average(log_n[ELL==0]),0.5*(perr_log_n[ELL==0]+merr_log_n[ELL==0]),mag_sph[ELL==0],
        	0.5*(merr_mag_sph[ELL==0] + perr_mag_sph[ELL==0]),log_n[ELL==0]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[ELL==0])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
	
	xx = np.arange(0.001,30,0.1)
	logxx = np.log10(xx)
	#print A[2], B[2]
	# calculates the prediction bands for the given input arrays
	lpb68,upb68,logxx = predband.predband(log_n[ELL==0]-np.average(log_n[ELL==0]),mag_sph[ELL==0],A[2],B[2],conf=0.68,x=logxx)
	#lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
	#lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
	yy = (A[2]*(logxx) + B[2])
	
	ax.plot(logxx+np.average(log_n[ELL==0]),yy, color='blue', linewidth=2.)
	# plots a shaded area containing the prediction band  
	plt.fill_between(logxx+np.average(log_n[ELL==0]), lpb68, upb68, alpha=0.15, facecolor='blue')
	ax.plot(logxx+np.average(log_n[ELL==0]),lpb68, color='blue')
	ax.plot(logxx+np.average(log_n[ELL==0]),upb68, color='blue')

	
        #ax.plot(logn,mag_sph, color='green', linewidth=2.)
        # plots a shaded area containing the prediction band  
        #plt.fill_between(logxx+np.average(logn), 10**lpb68, 10**upb68, alpha=0.15, facecolor='green')
        #ax.plot(logxx+np.average(logn),10**lpb68, color='green')
        #ax.plot(logxx+np.average(logn),10**upb68, color='green')
        #plt.fill_between(logxx+np.average(logn), 10**lpb95, 10**upb95, alpha=0.1, facecolor='green')
        #ax.plot(logxx+np.average(logn),10**lpb95, color='green')
        #ax.plot(logxx+np.average(logn),10**upb95, color='green')
        #plt.fill_between(logxx+np.average(logn), 10**lpb99, 10**upb99, alpha=0.05, facecolor='green')
        #ax.plot(logxx+np.average(logn),10**lpb99, color='green')
        #ax.plot(logxx+np.average(logn),10**upb99, color='green')
        plt.axis([np.log10(0.25),np.log10(18),-19.1,-27.9])
        plt.xlabel(r'$n^{\rm ' + axis + '}$', labelpad=20)
        plt.ylabel(r'$MAG_{\rm sph} \rm ~[mag]$', labelpad=20)
        plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
        #plt.savefig(path_scalrel_plots + 'mag_sph_vs_logn' + axis + '.pdf', format='pdf', dpi=1000)
        #plt.clf()
		
	
	

def log_mbh_vs_logn(axis):
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.ELLIPTICAL_my, anc.core, \
		pysres.log_n_' + axis + '_moffat_comb, \
		errV.perr_log_n, errV.merr_log_n \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	mbh = data[1].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[2].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[3].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
        ELL = data[4].astype(np.int)
	core = data[5].astype(np.int)
	log_n = data[6].astype(np.float)
	n = 10**log_n
	perr_log_n = data[7].astype(np.float)
	merr_log_n = data[8].astype(np.float)
		
	A,B,Aerr,Berr,covAB=bces.bces(log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),log_mbh,0.5*(merr_log_mbh + perr_log_mbh),log_n*[0.0])
	print '---------------------------------'
	print 'y = A*(x-<x>) + B '
	print '<x> =', np.average(log_n)
	print
	print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
	print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
	print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
	print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
	print '---------------------------------'
	
	xx = np.arange(0.001,30,0.1)
	logxx = np.log10(xx)
	#print A[2], B[2]
	# calculates the prediction bands for the given input arrays
	lpb68,upb68,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.68,x=logxx)
	lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
	lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
	yy = (A[2]*(logxx) + B[2])
	#print lpb, upb,xx
	
        fig, ax = plt.subplots()
        #ax.errorbar(log_n[ELL==1 ], mbh[ELL==1 ], xerr=[merr_log_n[ELL==1 ],perr_log_n[ELL==1 ]], 
	#	yerr=[merr_mbh[ELL==1 ],perr_mbh[ELL==1 ]], ecolor='gray', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        #ax.errorbar(log_n, mbh, xerr=[merr_log_n,perr_log_n], yerr=[merr_mbh,perr_mbh], ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.scatter(log_n[ELL==1] and (core==1)], mbh[(ELL==1) and (core==1)], c='black', s=60)
        ax.set_yscale('log')
        xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
       #ax.plot(logxx+np.average(log_n),10**yy, color='green', linewidth=2.)
       ## plots a shaded area containing the prediction band  
       #plt.fill_between(logxx+np.average(log_n), 10**lpb68, 10**upb68, alpha=0.15, facecolor='green')
       #ax.plot(logxx+np.average(log_n),10**lpb68, color='green')
       #ax.plot(logxx+np.average(log_n),10**upb68, color='green')
       #plt.fill_between(logxx+np.average(log_n), 10**lpb95, 10**upb95, alpha=0.1, facecolor='green')
       #ax.plot(logxx+np.average(log_n),10**lpb95, color='green')
       #ax.plot(logxx+np.average(log_n),10**upb95, color='green')
       #plt.fill_between(logxx+np.average(log_n), 10**lpb99, 10**upb99, alpha=0.05, facecolor='green')
       #ax.plot(logxx+np.average(log_n),10**lpb99, color='green')
       #ax.plot(logxx+np.average(log_n),10**upb99, color='green')
        plt.axis([np.log10(0.25),np.log10(18),10**5.5,10**11.2])
        plt.xlabel(r'$n^{\rm ' + axis + '}$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_log_n_' + axis + '.pdf', format='pdf', dpi=1000)
	plt.clf()
	
def log_mbh_vs_mag_sph():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	mbh = data[1].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[2].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[3].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	core = data[4].astype(np.int)
	mag_sph = data[5].astype(np.float)
	perr_mag_sph = data[6].astype(np.float)
	merr_mag_sph = data[7].astype(np.float)
		
        fig, ax = plt.subplots()

        print 'core-Sersic'
	A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==1]-np.average(mag_sph[core==1]),0.5*(perr_mag_sph[core==1]+merr_mag_sph[core==1]),
		log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),mag_sph[core==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[core==1])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-40,40,0.01)
        # calculates the prediction bands for the given input arrays
        lpb68,upb68,logxx = predband.predband(mag_sph[core==1]-np.average(mag_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.68,x=logxx)
        lpb95,upb95,logxx = predband.predband(mag_sph[core==1]-np.average(mag_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.95,x=logxx)
        lpb99,upb99,logxx = predband.predband(mag_sph[core==1]-np.average(mag_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
        #print lpb, upb,xx
        ax.plot(logxx+np.average(mag_sph[core==1]),10**yy, color='pink', linewidth=2.)
        # plots a shaded area containing the prediction band
        plt.fill_between(logxx+np.average(mag_sph[core==1]), 10**lpb68, 10**upb68, alpha=0.15, facecolor='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**lpb68, color='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**upb68, color='pink')
        plt.fill_between(logxx+np.average(mag_sph[core==1]), 10**lpb95, 10**upb95, alpha=0.1, facecolor='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**lpb95, color='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**upb95, color='pink')
        plt.fill_between(logxx+np.average(mag_sph[core==1]), 10**lpb99, 10**upb99, alpha=0.05, facecolor='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**lpb99, color='pink')
        ax.plot(logxx+np.average(mag_sph[core==1]),10**upb99, color='pink')
	
        print 'Sersic'
	A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==0]-np.average(mag_sph[core==0]),0.5*(perr_mag_sph[core==0]+merr_mag_sph[core==0]),
		log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),mag_sph[core==0]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[core==0])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-40,40,0.01)
        # calculates the prediction bands for the given input arrays
        lpb68,upb68,logxx = predband.predband(mag_sph[core==0]-np.average(mag_sph[core==0]),log_mbh[core==0],A[2],B[2],conf=0.68,x=logxx)
        lpb95,upb95,logxx = predband.predband(mag_sph[core==0]-np.average(mag_sph[core==0]),log_mbh[core==0],A[2],B[2],conf=0.95,x=logxx)
        lpb99,upb99,logxx = predband.predband(mag_sph[core==0]-np.average(mag_sph[core==0]),log_mbh[core==0],A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
        #print lpb, upb,xx
        ax.plot(logxx+np.average(mag_sph[core==0]),10**yy, color='cyan', linewidth=2.)
        # plots a shaded area containing the prediction band
        plt.fill_between(logxx+np.average(mag_sph[core==0]), 10**lpb68, 10**upb68, alpha=0.15, facecolor='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**lpb68, color='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**upb68, color='cyan')
        plt.fill_between(logxx+np.average(mag_sph[core==0]), 10**lpb95, 10**upb95, alpha=0.1, facecolor='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**lpb95, color='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**upb95, color='cyan')
        plt.fill_between(logxx+np.average(mag_sph[core==0]), 10**lpb99, 10**upb99, alpha=0.05, facecolor='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**lpb99, color='cyan')
        ax.plot(logxx+np.average(mag_sph[core==0]),10**upb99, color='cyan')
	
	
        ax.errorbar(mag_sph[core==1], mbh[core==1], xerr=[merr_mag_sph[core==1],perr_mag_sph[core==1]], yerr=[merr_mbh[core==1],perr_mbh[core==1]], ecolor='gray', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mag_sph[core==0], mbh[core==0], xerr=[merr_mag_sph[core==0],perr_mag_sph[core==0]], yerr=[merr_mbh[core==0],perr_mbh[core==0]], ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        #ax.scatter(mag_sph, mbh, c='black', s=60)
        ax.set_yscale('log')
        plt.axis([-19.01,-27.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph.pdf', format='pdf', dpi=1000)
	plt.clf()
	
def log_mbh_vs_mu_0(axis):
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, \
		pysres.log_n_' + axis + '_moffat_comb, pysres.mu_e_' + axis + '_moffat_comb, \
		errV.perr_mu_0, errV.merr_mu_0 \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	mbh = data[1].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[2].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[3].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	log_n = data[4].astype(np.float)
	n = 10**log_n
	mu_e = data[5].astype(np.float)
	perr_mu_0 = data[6].astype(np.float)
	merr_mu_0 = data[7].astype(np.float)
	
	# compute mu_0
	b = mu_e * [0.0]
	for i in range(len(b)):
		b[i] = b_n.computeb_n(n[i])
	
	mu_0 = mu_e - 2.5*b/np.log(10)
		
	A,B,Aerr,Berr,covAB=bces.bces(mu_0-np.average(mu_0),0.5*(merr_mu_0+perr_mu_0),log_mbh,0.5*(merr_log_mbh + perr_log_mbh),mu_0*[0.0])
	print '---------------------------------'
	print 'y = A*(x-<x>) + B '
	print '<x> =', np.average(mu_0)
	print
	print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
	print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
	print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
	print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
	print '---------------------------------'
	
	#xx = np.arange(0.001,30,0.1)
	logxx = np.arange(-40,40,0.1)
	#print A[2], B[2]
	# calculates the prediction bands for the given input arrays
	lpb68,upb68,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.68,x=logxx)
	lpb95,upb95,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.95,x=logxx)
	lpb99,upb99,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.99,x=logxx)
	yy = (A[2]*(logxx) + B[2])
	#print lpb, upb,xx
	
        fig, ax = plt.subplots()
        ax.errorbar(mu_0, mbh, xerr=[merr_mu_0,perr_mu_0], yerr=[merr_mbh,perr_mbh], ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        #plt.scatter(logn, mbh, c='black')
        ax.set_yscale('log')
        #xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        #xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        #ax.set_xticks(xticks)
        #ax.set_xticklabels(xticks_labels)
	ax.plot(logxx+np.average(mu_0),10**yy, color='green', linewidth=2.)
	# plots a shaded area containing the prediction band  
	plt.fill_between(logxx+np.average(mu_0), 10**lpb68, 10**upb68, alpha=0.15, facecolor='green')
	ax.plot(logxx+np.average(mu_0),10**lpb68, color='green')
	ax.plot(logxx+np.average(mu_0),10**upb68, color='green')
	plt.fill_between(logxx+np.average(mu_0), 10**lpb95, 10**upb95, alpha=0.1, facecolor='green')
	ax.plot(logxx+np.average(mu_0),10**lpb95, color='green')
	ax.plot(logxx+np.average(mu_0),10**upb95, color='green')
	plt.fill_between(logxx+np.average(mu_0), 10**lpb99, 10**upb99, alpha=0.05, facecolor='green')
	ax.plot(logxx+np.average(mu_0),10**lpb99, color='green')
	ax.plot(logxx+np.average(mu_0),10**upb99, color='green')
        plt.axis([-10.1,-34.9,10**5.5,10**11.2])
        plt.xlabel(r'$\mu_0^{\rm ' + axis + r'} \rm ~[mag~arcsec^{-2}]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.17)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mu_0_' + axis + '.pdf', format='pdf', dpi=1000)
	#plt.clf()
	
		

#def mbh_vs_log_r_e():
       #connection = sql3.connect(dbname)
       #cur = connection.cursor()
       #
       #getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, \
       #	pysres.log_r_e_eq_moffat_comb, \
       #	errC.log_r_e_std, errC.log_r_e_num_measurements, \
       #	errV.err_n_vote \
       #	FROM Ancillary AS anc \
       #	JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
       #	JOIN ErrorsComparison as errC ON anc.gal_id = errC.gal_id \
       #	JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
       #	WHERE anc.fit1D_done = 1;'
       #cur.execute(getdata_query)
       #datalist = cur.fetchall()
       #data= np.asarray(datalist).transpose()
       ##print data
       #mbh = data[1].astype(np.float)
       #logmbh = np.log10(mbh)
       #perr_mbh = data[2].astype(np.float)
       #perr_logmbh = np.log10(1 + perr_mbh/mbh)
       #merr_mbh = data[3].astype(np.float)
       #merr_logmbh = -np.log10(1 - merr_mbh/mbh)
       #log_r_e = data[4].astype(np.float)
       #r_e = 10**log_r_e
       #if (error_type=='comparison'):
       #	err_log_r_e = data[5].astype(np.float)
       #	num_measurements_log_r_e = data[6].astype(np.int)
       #
       #	err_log_r_e_without0 = err_log_r_e
       #	for i in range(len(err_log_r_e)):
       #		if num_measurements_log_r_e[i] == 1:
       #			np.delete(err_log_r_e_without0, i) 
       #
       #	# compute average error
       #	#average_err_logn = np.average(err_logn_without0)
       #	#average_err_n_percentage = (10**average_err_logn - 1)
       #	#print len(num_measurements_logn[num_measurements_logn>1]), average_err_logn, average_err_n_percentage
       #	# compute median error
       #	median_err_log_r_e = np.median(err_log_r_e_without0)
       #	median_err_r_e_percentage = (10**median_err_log_r_e - 1)
       #	print 'median error on log(R_e) = ', median_err_log_r_e
       #	print 'median percentage error on R_e = ', median_err_r_e_percentage
       #	# assign default error to measurements with no counterpart in literature
       #	for i in range(len(err_log_r_e)):
       #		if num_measurements_log_r_e[i] == 1:
       #			err_log_r_e[i] = median_err_log_r_e
       #	perr_r_e = (10**err_log_r_e - 1)*r_e
       #	merr_r_e = (10**(-err_log_r_e) + 1)*(-r_e)
       #
       #	histo_err_log_r_e(err_log_r_e_without0,median_err_log_r_e)
       #elif error_type=='vote':
       #	err_n_vote = data[7].astype(np.int)
       #	err_log_r_e = log_r_e*0.0 + np.log10(1 + default_err_r_e_percentage + 0.15*(err_n_vote-1))		
       #
       #A,B,Aerr,Berr,covAB=bces.bces(log_r_e-np.average(log_r_e),err_log_r_e,logmbh,0.5*(merr_logmbh + perr_logmbh),log_r_e*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_r_e)
       #print
       #print 'OLS(Y|X)    A =', A[0], '+-', Aerr[0], '   B = ', B[0], '+-', Berr[0]
       #print 'OLS(X|Y)    A =', A[1], '+-', Aerr[1], '   B = ', B[1], '+-', Berr[1]
       #print 'bisector    A =', A[2], '+-', Aerr[2], '   B = ', B[2], '+-', Berr[2]
       #print 'orthogonal  A =', A[3], '+-', Aerr[3], '   B = ', B[3], '+-', Berr[3]
       #print '---------------------------------'
       #
       #xx = np.arange(0.001,100,0.1)
       #logxx = np.log10(xx)
       ##print A[2], B[2]
       ## calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(log_r_e-np.average(log_r_e),logmbh,A[2],B[2],conf=0.68,x=logxx)
       #lpb95,upb95,logxx = predband.predband(log_r_e-np.average(log_r_e),logmbh,A[2],B[2],conf=0.95,x=logxx)
       #lpb99,upb99,logxx = predband.predband(log_r_e-np.average(log_r_e),logmbh,A[2],B[2],conf=0.99,x=logxx)
       #yy = (A[2]*(logxx) + B[2])
       ##print lpb, upb,xx
       #
       #
       #fig, ax = plt.subplots()
       #ax.errorbar(log_r_e, mbh, xerr=err_log_r_e, yerr=[merr_mbh,perr_mbh], ecolor='black', fmt=None, elinewidth=1.5, capthick=1.5) 
       ##plt.scatter(logn, mbh, c='black')
       #ax.set_yscale('log')
       #xticks = np.log10(np.asarray([0.1, 0.5, 1., 2., 3., 10., 20., 50.]))
       #xticks_labels = ['$0.1$','$0.5$','$1$','$2$','$3$', '$10$','$20$','$50$']
       #ax.set_xticks(xticks)
       #ax.set_xticklabels(xticks_labels)
       #ax.plot(logxx+np.average(log_r_e),10**yy, color='green', linewidth=2.)
       ## plots a shaded area containing the prediction band  
       #plt.fill_between(logxx+np.average(log_r_e), 10**lpb68, 10**upb68, alpha=0.15, facecolor='green')
       #ax.plot(logxx+np.average(log_r_e),10**lpb68, color='green')
       #ax.plot(logxx+np.average(log_r_e),10**upb68, color='green')
       #plt.fill_between(logxx+np.average(log_r_e), 10**lpb95, 10**upb95, alpha=0.1, facecolor='green')
       #ax.plot(logxx+np.average(log_r_e),10**lpb95, color='green')
       #ax.plot(logxx+np.average(log_r_e),10**upb95, color='green')
       #plt.fill_between(logxx+np.average(log_r_e), 10**lpb99, 10**upb99, alpha=0.05, facecolor='green')
       #ax.plot(logxx+np.average(log_r_e),10**lpb99, color='green')
       #ax.plot(logxx+np.average(log_r_e),10**upb99, color='green')
       #plt.axis([np.log10(0.035),np.log10(99),10**5.5,10**11.2])
       #plt.xlabel(r'$R_{\rm e}^{\rm eq} \rm ~ [kpc]$', labelpad=20)
       #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
       #plt.subplots_adjust(left=0.15,bottom=0.15)
       ##plt.show()
       #plt.savefig(path_scalrel_plots + 'mbh_vs_log_r_e_' + error_type + '.pdf', format='pdf', dpi=1000)
       #plt.clf()
		


def main():
	log_mbh_vs_logn('maj')	
	#log_mbh_vs_logn('eq')		
	#log_mbh_vs_mu_0('maj')
	#log_mbh_vs_mu_0('eq')
	#log_mbh_vs_mag_sph()
	#mbh_vs_logn('eq')		
	#mbh_vs_logr_e('comparison')
	#mbh_vs_logr_e()
	#mag_sph_vs_logn('maj')
	#mag_sph_vs_logn('eq')
main()
