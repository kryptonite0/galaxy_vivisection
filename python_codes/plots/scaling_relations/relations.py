import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy
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

def mag_tot_GS13_vs_me():
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.KMAG_tot, \
		pysres.mag_tot_eq_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		WHERE anc.fit1D_done = 1 \
		ORDER BY pysres.mag_tot_eq_moffat_comb;'
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	kmag_tot = data[1].astype(np.float)
	mag_tot = data[2].astype(np.float)
	
	for a,b,c in zip(gal_id,mag_tot,kmag_tot):	
		print a,b,c
	
	fig, ax = plt.subplots()
	ax.scatter(mag_tot,kmag_tot)
	ax.plot(np.arange(-30,-19,1),np.arange(-30,-19,1)+0.27, ls='--')
	ax.axis([-21.01,-27.99,-21.01,-27.99])
	plt.xlabel(r'$3.6{\rm~\mu m}~MAG_{\rm tot} \rm ~[mag]$', labelpad=20)
        plt.ylabel(r'$K$-band $MAG_{\rm tot} \rm ~[mag]$', labelpad=20)
        plt.subplots_adjust(left=0.20,bottom=0.20)
	#plt.show()
	plt.savefig(path_scalrel_plots + 'mag_tot_GS13_vs_me.pdf', format='pdf', dpi=1000)

def BTratios_GS13_vs_me():
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.morphtype, \
		pysres.mag_sph_eq_moffat_comb, pysres.mag_tot_eq_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		WHERE anc.fit1D_done = 1 \
		ORDER BY pysres.mag_sph_eq_moffat_comb;'
		
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	morphtype = data[1]
	mag_sph = data[2].astype(np.float)	
	mag_tot = data[3].astype(np.float)
	
	BT_my = 10**(-0.4*(mag_sph))/10**(-0.4*(mag_tot))
	
	BT_GS13 = mag_sph*[0.0] + 1.0
	BT_GS13[morphtype=='S0'] = 0.49
	BT_GS13[morphtype=='SB0'] = 0.49
	BT_GS13[morphtype=='S0a'] = 0.49
	BT_GS13[morphtype=='SB0a'] = 0.49
	BT_GS13[morphtype=='Sa'] = 0.46
	BT_GS13[morphtype=='SBa'] = 0.46
	BT_GS13[morphtype=='Sab'] = 0.29
	BT_GS13[morphtype=='SBab'] = 0.29
	BT_GS13[morphtype=='Sb'] = 0.25
	BT_GS13[morphtype=='SBb'] = 0.25
	BT_GS13[morphtype=='Sbc'] = 0.15
	BT_GS13[morphtype=='SBbc'] = 0.15		
	BT_GS13[morphtype=='Sc'] = 0.09
	BT_GS13[morphtype=='SBc'] = 0.09
	BT_GS13[morphtype=='Scd'] = 0.06
	BT_GS13[morphtype=='SBcd'] = 0.06
	
	for a,b,c,d,e in zip(gal_id,mag_sph,morphtype,BT_my,BT_GS13):
		print a,round(b,2),c,round(d,2),e 
		
	diff = BT_my - BT_GS13	
	
	fig, ax = plt.subplots()
	ax.scatter(mag_sph,diff)
	plt.axis([-19.01,-27.99,-1,1])
	plt.xlabel(r'$3.6{\rm~\mu m}~MAG_{\rm sph} \rm ~[mag]$', labelpad=20)
        plt.ylabel(r'$B/T_{my} - B/T_{GS13}$', labelpad=20)
        plt.subplots_adjust(left=0.20,bottom=0.20)
	#plt.show()
	plt.savefig(path_scalrel_plots + 'BTratios_GS13_vs_me.pdf', format='pdf', dpi=1000)
	
	
def sani_vs_me():
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id,  \
		pysres.mag_sph_eq_moffat_comb, \
		pyssani.mag_sph \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN LiteratureDecompositionsSanietal2011PhysicalUnits AS pyssani ON anc.gal_id = pyssani.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	mag_sph_my = data[1].astype(np.float)	
	mag_sph_sani = data[2].astype(np.float)	
	
	fig, ax = plt.subplots()
	ax.scatter(mag_sph_my, mag_sph_sani)
	ax.plot(np.arange(-28,-19,1),np.arange(-28,-19,1))
	plt.axis([-19.01,-27.99,-19.01,-27.99])
	plt.show()
		

def KMAG_sph_vs_mag_sph():
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.KMAG_sph, anc.ELLIPTICAL_my, anc.core, \
		pysres.mag_sph_eq_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	KMAG_sph = data[1].astype(np.float)
        ELL = data[2].astype(np.int)
	core = data[3].astype(np.int)
	mag_sph = data[4].astype(np.float)
	
	fig, ax = plt.subplots()
	ax.scatter(mag_sph[ELL==1], KMAG_sph[ELL==1], color='red', s=50)	
	ax.scatter(mag_sph[ELL==0], KMAG_sph[ELL==0], color='blue', s=50)
	ax.scatter(mag_sph[core==1], KMAG_sph[core==1], color='white', s=20)	
	
        full_names = {'n4388' : 'NGC 4388', 'n3585' : 'NGC 3585', 'n3607' : 'NGC 3607', 'n3414' : 'NGC 3414', 'n3115' : 'NGC 3115', 'n3377' : 'NGC 3377'}
	
	for g, m, k in zip(gal_id, mag_sph, KMAG_sph):
		if g in ['n4388']:
			ax.annotate(full_names[g], xy = (m, k), xytext = (0, 20),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'blue', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
		elif g in ['n3607']:	
        		ax.annotate(full_names[g], xy = (m, k), xytext = (110, -30),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'red', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
		elif g in ['n3414']:	
        		ax.annotate(full_names[g], xy = (m, k), xytext = (120, -35),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'red', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
		elif g in ['n3585']:	
        		ax.annotate(full_names[g], xy = (m, k), xytext = (100, 10),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'red', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
		elif g in ['n3115']:	
        		ax.annotate(full_names[g], xy = (m, k), xytext = (80, -80),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'red', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
		elif g in ['n3377']:	
        		ax.annotate(full_names[g], xy = (m, k), xytext = (40, -105),
        			textcoords = 'offset points', ha = 'right', va = 'bottom',
        			bbox = dict(boxstyle = 'round,pad=0.2', fc = 'red', alpha = 0.1),
        			arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))
				
	xx = np.arange(-50,0,1)
	ax.plot(xx,xx+0.27, color='k', ls='--')

        plt.axis([-19.01,-27.99,-19.01,-27.99])
        plt.xlabel(r'$3.6{\rm~\mu m}~MAG_{\rm sph} \rm ~[mag]$', labelpad=20)
        plt.ylabel(r'$K-$band $MAG_{\rm sph} \rm ~[mag]$', labelpad=20)
        plt.subplots_adjust(left=0.15,bottom=0.15)
	#plt.show()
	plt.savefig(path_scalrel_plots + 'KMAG_sph_vs_mag_sph.pdf', format='pdf', dpi=1000)

def log_mbh_vs_KMAG_sph():
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.KMAG_sph, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, anc.ELLIPTICAL_my \
		FROM Ancillary AS anc \
		WHERE anc.source = "GS13" AND anc.KMAG_sph>-99;'
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	KMAG_sph = data[1].astype(np.float)
        mbh = data[2].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[3].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[4].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	core = data[5].astype(np.int)
        ELL = data[6].astype(np.int)
	
	err_KMAG_sph = KMAG_sph*[0.0] + 0.25
	err_KMAG_sph[ELL==0] = 0.75
	
	fig, ax = plt.subplots()
        
	### fit using FITEXY ###
        print 'core-Sersic'
        A,B = fitexy.bisect_modfitexy(KMAG_sph[core==1]-np.average(KMAG_sph[core==1]), err_KMAG_sph[core==1],
        	log_mbh[core==1], 0.5*(perr_log_mbh[core==1] + merr_log_mbh[core==1]))
        logxx = np.arange(-50,50,0.1)
        # plot y|x relation
        y_1 = A[0] + B[0]*logxx
       # ax.plot(logxx+np.average(KMAG_sph[core==1]), 10**y_1, ls='--', color='gray', linewidth=2.)
        # plot x|y relation
        y_2 = A[1] + B[1]*logxx
       # ax.plot(logxx+np.average(KMAG_sph[core==1]), 10**y_2, ls='--', color='gray', linewidth=2.)
        # plot bisectore relation
        y_bisec = A[2] + B[2]*logxx
        ax.plot(logxx+np.average(KMAG_sph[core==1]), 10**y_bisec, ls='-', color='gray', linewidth=2.)
	slope_cS = B[2]
      
        print 'Sersic'
        A,B = fitexy.bisect_modfitexy(KMAG_sph[core==0]-np.average(KMAG_sph[core==0]), err_KMAG_sph[core==0],
        	log_mbh[core==0], 0.5*(perr_log_mbh[core==0] + merr_log_mbh[core==0]))
        logxx = np.arange(-50,50,0.1)
        # plot y|x relation
        y_1 = A[0] + B[0]*logxx
       # ax.plot(logxx+np.average(KMAG_sph[core==0]), 10**y_1, ls='--', color='black', linewidth=2.)
        # plot x|y relation
        y_2 = A[1] + B[1]*logxx
       # ax.plot(logxx+np.average(KMAG_sph[core==0]), 10**y_2, ls='--', color='black', linewidth=2.)
        # plot bisectore relation
        y_bisec = A[2] + B[2]*logxx
        ax.plot(logxx+np.average(KMAG_sph[core==0]), 10**y_bisec, ls='-', color='black', linewidth=2.)
	slope_S = B[2]
	print 'slope_cS=', slope_cS
	print 'slope_S=', slope_S
	
	### fit using BCES ####
        print 'core-Sersic'
        A,B,Aerr,Berr,covAB=bces.bces(KMAG_sph[core==1]-np.average(KMAG_sph[core==1]), err_KMAG_sph[core==1], log_mbh[core==1],
        	0.5*(perr_log_mbh[core==1] + merr_log_mbh[core==1]),KMAG_sph[core==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(KMAG_sph[core==1])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-50,50,0.1)
	y_bcesbisec = (A[2]*(logxx) + B[2])
	ax.plot(logxx+np.average(KMAG_sph[core==1]),10**y_bcesbisec, color='red', linewidth=2.)
	
        print 'Sersic'
        A,B,Aerr,Berr,covAB=bces.bces(KMAG_sph[core==0]-np.average(KMAG_sph[core==0]), err_KMAG_sph[core==0], log_mbh[core==0],
        	0.5*(perr_log_mbh[core==0] + merr_log_mbh[core==0]),KMAG_sph[core==0]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(KMAG_sph[core==0])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-50,50,0.1)
	y_bcesbisec = (A[2]*(logxx) + B[2])
	ax.plot(logxx+np.average(KMAG_sph[core==0]),10**y_bcesbisec, color='blue', linewidth=2.)
	


	ax.set_yscale('log')
        ax.errorbar(KMAG_sph[core==1], mbh[core==1], xerr=err_KMAG_sph[core==1], yerr=[merr_mbh[core==1],perr_mbh[core==1]], ecolor='gray', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(KMAG_sph[core==0], mbh[core==0], xerr=err_KMAG_sph[core==0], yerr=[merr_mbh[core==0],perr_mbh[core==0]], ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
	plt.axis([-18.01,-27.99,10**5.001,10**10.999])
	plt.show()


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
	
	### fit using BCES ####
       #print 'Ellipticals'
       ##print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
       #A,B,Aerr,Berr,covAB=bces.bces(log_n[ELL==1]-np.average(log_n[ELL==1]),0.5*(perr_log_n[ELL==1]+merr_log_n[ELL==1]),mag_sph[ELL==1],
       #	0.5*(merr_mag_sph[ELL==1] + perr_mag_sph[ELL==1]),log_n[ELL==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n[ELL==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #xx = np.arange(0.001,30,0.1)
       #logxx = np.log10(xx)
       ##print A[2], B[2]
       ## calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(log_n[ELL==1]-np.average(log_n[ELL==1]),mag_sph[ELL==1],A[2],B[2],conf=0.68,x=logxx)
       ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       ##lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       #yy = (A[2]*(logxx) + B[2])
       #
       #ax.plot(logxx+np.average(log_n[ELL==1]),yy, color='red', linewidth=2.)
       ## plots a shaded area containing the prediction band  
       #plt.fill_between(logxx+np.average(log_n[ELL==1]), lpb68, upb68, alpha=0.15, facecolor='red')
       #ax.plot(logxx+np.average(log_n[ELL==1]),lpb68, color='red')
       #ax.plot(logxx+np.average(log_n[ELL==1]),upb68, color='red')
       #
       #print 'Bulges'
       #A,B,Aerr,Berr,covAB=bces.bces(log_n[ELL==0]-np.average(log_n[ELL==0]),0.5*(perr_log_n[ELL==0]+merr_log_n[ELL==0]),mag_sph[ELL==0],
       #	0.5*(merr_mag_sph[ELL==0] + perr_mag_sph[ELL==0]),log_n[ELL==0]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n[ELL==0])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #xx = np.arange(0.001,30,0.1)
       #logxx = np.log10(xx)
       ##print A[2], B[2]
       ## calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(log_n[ELL==0]-np.average(log_n[ELL==0]),mag_sph[ELL==0],A[2],B[2],conf=0.68,x=logxx)
       ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       ##lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       #yy = (A[2]*(logxx) + B[2])
       #
       #ax.plot(logxx+np.average(log_n[ELL==0]),yy, color='blue', linewidth=2.)
       ## plots a shaded area containing the prediction band  
       #plt.fill_between(logxx+np.average(log_n[ELL==0]), lpb68, upb68, alpha=0.15, facecolor='blue')
       #ax.plot(logxx+np.average(log_n[ELL==0]),lpb68, color='blue')
       #ax.plot(logxx+np.average(log_n[ELL==0]),upb68, color='blue')
       
        ### fit using FITEXY ###
        print 'Ellipticals'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_n[ELL==1]-np.average(log_n[ELL==1]), 0.5*(perr_log_n[ELL==1]+merr_log_n[ELL==1]),
        	mag_sph[ELL==1], 0.5*(merr_mag_sph[ELL==1] + perr_mag_sph[ELL==1]))
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        # plot bisector relation
        y_bisec = A + B*logxx
        ax.plot(logxx+np.average(log_n[ELL==1]), y_bisec, ls='-', color='red', linewidth=2.)
	
        ### fit using FITEXY ###
	print 'Bulges'
	A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_n[ELL==0]-np.average(log_n[ELL==0]), 0.5*(perr_log_n[ELL==0]+merr_log_n[ELL==0]),
		mag_sph[ELL==0], 0.5*(merr_mag_sph[ELL==0] + perr_mag_sph[ELL==0]))
	xx = np.arange(0.001,30,0.1)
	logxx = np.log10(xx)
	# plot bisector relation
	y_bisec = A + B*logxx
	ax.plot(logxx+np.average(log_n[ELL==0]), y_bisec, ls='-', color='blue', linewidth=2.)
	
	#print A,B
	
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
        #plt.show()
        plt.savefig(path_scalrel_plots + 'mag_sph_vs_logn' + axis + '.pdf', format='pdf', dpi=1000)
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
        #ax.scatter(log_n[ELL==1] and (core==1)], mbh[(ELL==1) and (core==1)], c='black', s=60)
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
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my \
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
	ELL = data[8].astype(np.int)
	
	ELL_core = ELL*[0]
	ELL_sersic = ELL*[0]
	BUL_core = ELL*[0]
	BUL_sersic = ELL*[0]
	
	for i in range(len(ELL)):
		if ELL[i]==1 and core[i]==1:
			ELL_core[i] = 1
		elif ELL[i]==1 and core[i]==0:	
			ELL_sersic[i] = 1
		elif ELL[i]==0 and core[i]==1:	
			BUL_core[i] = 1
		elif ELL[i]==0 and core[i]==0:	
			BUL_sersic[i] = 1	
	
        fig, ax = plt.subplots()

      ##print 'ELL'
      ##A,B,Aerr,Berr,covAB=bces.bces(mag_sph[ELL==1]-np.average(mag_sph[ELL==1]),0.5*(perr_mag_sph[ELL==1]+merr_mag_sph[ELL==1]),
      ##	log_mbh[ELL==1],0.5*(merr_log_mbh[ELL==1] + perr_log_mbh[ELL==1]),mag_sph[ELL==1]*[0.0])
      ##print '---------------------------------'
      ##print 'y = A*(x-<x>) + B '
      ##print '<x> =', np.average(mag_sph[ELL==1])
      ##print
      ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      ##print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      ##print '---------------------------------'
      ##
      ##logxx = np.arange(-40,40,0.01)
      ##yy = (A[2]*(logxx) + B[2])
      ##ax.plot(logxx+np.average(mag_sph[ELL==1]),10**yy, color='red', linewidth=2.)
      ##
      ##print 'BUL'
      ##A,B,Aerr,Berr,covAB=bces.bces(mag_sph[ELL==0]-np.average(mag_sph[ELL==0]),0.5*(perr_mag_sph[ELL==0]+merr_mag_sph[ELL==0]),
      ##	log_mbh[ELL==0],0.5*(merr_log_mbh[ELL==0] + perr_log_mbh[ELL==0]),mag_sph[ELL==0]*[0.0])
      ##print '---------------------------------'
      ##print 'y = A*(x-<x>) + B '
      ##print '<x> =', np.average(mag_sph[ELL==0])
      ##print
      ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      ##print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      ##print '---------------------------------'
      ##
      ##yy = (A[2]*(logxx) + B[2])
      ##ax.plot(logxx+np.average(mag_sph[ELL==0]),10**yy, color='blue', linewidth=2.)
	
        ### fit using FITEXY ###
        print 'core-Sersic'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==1]-np.average(mag_sph[core==1]), 0.5*(perr_mag_sph[core==1]+merr_mag_sph[core==1]),
        	log_mbh[core==1], 0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
        logxx = np.arange(-40,40,0.01)
        # plot bisector relation
        y_bisec = A + B*logxx
        ax.plot(logxx+np.average(mag_sph[core==1]), 10**y_bisec, ls='--', color='k', linewidth=2.)
        label_cS = r'$B_{\rm c-Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
        ax.text(-19.5, 10**10.5, label_cS, fontsize=20)
      
        ### fit using FITEXY ###
        print 'Sersic'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==0]-np.average(mag_sph[core==0]), 0.5*(perr_mag_sph[core==0]+merr_mag_sph[core==0]),
        	log_mbh[core==0], 0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]))
        logxx = np.arange(-40,40,0.01)
        # plot bisector relation
        y_bisec = A + B*logxx
        ax.plot(logxx+np.average(mag_sph[core==0]), 10**y_bisec, ls='-', color='k', linewidth=2.)
        label_S = r'$B_{\rm Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
        ax.text(-19.5, 10**10, label_S, fontsize=20)
	
        ### fit using FITEXY ###
        print 'ELL'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[ELL==1]-np.average(mag_sph[ELL==1]), 0.5*(perr_mag_sph[ELL==1]+merr_mag_sph[ELL==1]),
        	log_mbh[ELL==1], 0.5*(merr_log_mbh[ELL==1] + perr_log_mbh[ELL==1]))
        logxx = np.arange(-40,40,0.01)
        # plot bisector relation
        y_bisec = A + B*logxx
        ax.plot(logxx+np.average(mag_sph[ELL==1]), 10**y_bisec, ls='-', color='red', linewidth=2.)
        label_ell = r'$B_{\rm ell} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
        ax.text(-24, 10**6.7, label_ell, fontsize=20)
      
        ### fit using FITEXY ###
        print 'BUL'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[ELL==0]-np.average(mag_sph[ELL==0]), 0.5*(perr_mag_sph[ELL==0]+merr_mag_sph[ELL==0]),
        	log_mbh[ELL==0], 0.5*(merr_log_mbh[ELL==0] + perr_log_mbh[ELL==0]))
        logxx = np.arange(-40,40,0.01)
        # plot bisector relation
        y_bisec = A + B*logxx
        ax.plot(logxx+np.average(mag_sph[ELL==0]), 10**y_bisec, ls='-', color='blue', linewidth=2.)
        label_bul = r'$B_{\rm bul} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
        ax.text(-24, 10**6, label_bul, fontsize=20)
	
        ax.errorbar(mag_sph[ELL_core==1], mbh[ELL_core==1], xerr=[merr_mag_sph[ELL_core==1],perr_mag_sph[ELL_core==1]], yerr=[merr_mbh[ELL_core==1],perr_mbh[ELL_core==1]], ecolor='red', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mag_sph[ELL_sersic==1], mbh[ELL_sersic==1], xerr=[merr_mag_sph[ELL_sersic==1],perr_mag_sph[ELL_sersic==1]], yerr=[merr_mbh[ELL_sersic==1],perr_mbh[ELL_sersic==1]], ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mag_sph[BUL_core==1], mbh[BUL_core==1], xerr=[merr_mag_sph[BUL_core==1],perr_mag_sph[BUL_core==1]], yerr=[merr_mbh[BUL_core==1],perr_mbh[BUL_core==1]], ecolor='blue', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mag_sph[BUL_sersic==1], mbh[BUL_sersic==1], xerr=[merr_mag_sph[BUL_sersic==1],perr_mag_sph[BUL_sersic==1]], yerr=[merr_mbh[BUL_sersic==1],perr_mbh[BUL_sersic==1]], ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.set_yscale('log')
        plt.axis([-19.01,-27.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph.pdf', format='pdf', dpi=1000)
	
def log_mbh_vs_mag_tot():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, pysres.mag_tot_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my \
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
	mag_tot = data[6].astype(np.float)
	perr_mag_sph = data[7].astype(np.float)
	merr_mag_sph = data[8].astype(np.float)
	ELL = data[9].astype(np.int)
	
	ELL_core = ELL*[0]
	ELL_sersic = ELL*[0]
	BUL_core = ELL*[0]
	BUL_sersic = ELL*[0]
	
	for i in range(len(ELL)):
		if ELL[i]==1 and core[i]==1:
			ELL_core[i] = 1
		elif ELL[i]==1 and core[i]==0:	
			ELL_sersic[i] = 1
		elif ELL[i]==0 and core[i]==1:	
			BUL_core[i] = 1
		elif ELL[i]==0 and core[i]==0:	
			BUL_sersic[i] = 1	
	
        fig, ax = plt.subplots()

      ##print 'ELL'
      ##A,B,Aerr,Berr,covAB=bces.bces(mag_sph[ELL==1]-np.average(mag_sph[ELL==1]),0.5*(perr_mag_sph[ELL==1]+merr_mag_sph[ELL==1]),
      ##	log_mbh[ELL==1],0.5*(merr_log_mbh[ELL==1] + perr_log_mbh[ELL==1]),mag_sph[ELL==1]*[0.0])
      ##print '---------------------------------'
      ##print 'y = A*(x-<x>) + B '
      ##print '<x> =', np.average(mag_sph[ELL==1])
      ##print
      ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      ##print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      ##print '---------------------------------'
      ##
      ##logxx = np.arange(-40,40,0.01)
      ##yy = (A[2]*(logxx) + B[2])
      ##ax.plot(logxx+np.average(mag_sph[ELL==1]),10**yy, color='red', linewidth=2.)
      ##
      ##print 'BUL'
      ##A,B,Aerr,Berr,covAB=bces.bces(mag_sph[ELL==0]-np.average(mag_sph[ELL==0]),0.5*(perr_mag_sph[ELL==0]+merr_mag_sph[ELL==0]),
      ##	log_mbh[ELL==0],0.5*(merr_log_mbh[ELL==0] + perr_log_mbh[ELL==0]),mag_sph[ELL==0]*[0.0])
      ##print '---------------------------------'
      ##print 'y = A*(x-<x>) + B '
      ##print '<x> =', np.average(mag_sph[ELL==0])
      ##print
      ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      ##print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      ##print '---------------------------------'
      ##
      ##yy = (A[2]*(logxx) + B[2])
      ##ax.plot(logxx+np.average(mag_sph[ELL==0]),10**yy, color='blue', linewidth=2.)
	
       #### fit using FITEXY ###
       #print 'core-Sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==1]-np.average(mag_sph[core==1]), 0.5*(perr_mag_sph[core==1]+merr_mag_sph[core==1]),
       #	log_mbh[core==1], 0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[core==1]), 10**y_bisec, ls='--', color='k', linewidth=2.)
       #label_cS = r'$B_{\rm c-Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       #ax.text(-19.5, 10**10.5, label_cS, fontsize=20)
       #
       #### fit using FITEXY ###
       #print 'Sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==0]-np.average(mag_sph[core==0]), 0.5*(perr_mag_sph[core==0]+merr_mag_sph[core==0]),
       #	log_mbh[core==0], 0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[core==0]), 10**y_bisec, ls='-', color='k', linewidth=2.)
       #label_S = r'$B_{\rm Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       #ax.text(-19.5, 10**10, label_S, fontsize=20)
       #
       #### fit using FITEXY ###
       #print 'ELL'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[ELL==1]-np.average(mag_sph[ELL==1]), 0.5*(perr_mag_sph[ELL==1]+merr_mag_sph[ELL==1]),
       #	log_mbh[ELL==1], 0.5*(merr_log_mbh[ELL==1] + perr_log_mbh[ELL==1]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[ELL==1]), 10**y_bisec, ls='-', color='red', linewidth=2.)
       #label_ell = r'$B_{\rm ell} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       #ax.text(-24, 10**6.7, label_ell, fontsize=20)
       #
       #### fit using FITEXY ###
       #print 'BUL'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[ELL==0]-np.average(mag_sph[ELL==0]), 0.5*(perr_mag_sph[ELL==0]+merr_mag_sph[ELL==0]),
       #	log_mbh[ELL==0], 0.5*(merr_log_mbh[ELL==0] + perr_log_mbh[ELL==0]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[ELL==0]), 10**y_bisec, ls='-', color='blue', linewidth=2.)
       #label_bul = r'$B_{\rm bul} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       #ax.text(-24, 10**6, label_bul, fontsize=20)
       #
       #ax.errorbar(mag_sph[ELL_core==1], mbh[ELL_core==1], xerr=[merr_mag_sph[ELL_core==1],perr_mag_sph[ELL_core==1]], yerr=[merr_mbh[ELL_core==1],perr_mbh[ELL_core==1]], ecolor='red', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[ELL_sersic==1], mbh[ELL_sersic==1], xerr=[merr_mag_sph[ELL_sersic==1],perr_mag_sph[ELL_sersic==1]], yerr=[merr_mbh[ELL_sersic==1],perr_mbh[ELL_sersic==1]], ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[BUL_core==1], mbh[BUL_core==1], xerr=[merr_mag_sph[BUL_core==1],perr_mag_sph[BUL_core==1]], yerr=[merr_mbh[BUL_core==1],perr_mbh[BUL_core==1]], ecolor='blue', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[BUL_sersic==1], mbh[BUL_sersic==1], xerr=[merr_mag_sph[BUL_sersic==1],perr_mag_sph[BUL_sersic==1]], yerr=[merr_mbh[BUL_sersic==1],perr_mbh[BUL_sersic==1]], ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
        ax.scatter(mag_tot[core==1],mbh[core==1], marker='*', color='k')
        ax.scatter(mag_tot[core==0],mbh[core==0], marker='o', color='k')
        ax.scatter(mag_sph[core==1],mbh[core==1], marker='*', color='red')
        ax.scatter(mag_sph[core==0],mbh[core==0], marker='o', color='red')
	ax.set_yscale('log')
        plt.axis([-19.01,-28.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm tot}\rm~[mag]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_tot.pdf', format='pdf', dpi=1000)
	
def mag_tot_vs_sigma():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.core, anc.ELLIPTICAL_my, anc.sigma, \
		pysres.mag_sph_eq_moffat_comb, pysres.mag_tot_eq_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	core = data[1].astype(np.int)
	ELL = data[2].astype(np.int)
	sigma = data[3].astype(np.float)
	mag_sph = data[4].astype(np.float)
	mag_tot = data[5].astype(np.float)
	
	
        fig, ax = plt.subplots()

        ax.scatter(sigma[ELL==1], mag_tot[ELL==1], marker='o', color='red', s=30)
        ax.scatter(sigma[ELL==0], mag_tot[ELL==0], marker='o', color='blue', s=30)
	ax.set_xscale('log')
	
        plt.axis([70,400,-19.01,-28.99])
        plt.xlabel(r'$\sigma$', labelpad=20)
        plt.ylabel(r'$magtot$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_tot.pdf', format='pdf', dpi=1000)
	
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
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mu_0_' + axis + '.pdf', format='pdf', dpi=1000)
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
	#log_mbh_vs_logn('maj')	
	#log_mbh_vs_logn('eq')		
	log_mbh_vs_mu_0('maj')
	#log_mbh_vs_mu_0('eq')
	#log_mbh_vs_mag_sph()
	#log_mbh_vs_mag_tot()
	#mbh_vs_logn('eq')		
	#mbh_vs_logr_e('comparison')
	#mbh_vs_logr_e()
	#mag_sph_vs_logn('maj')
	#mag_sph_vs_logn('eq')
	#KMAG_sph_vs_mag_sph()
	#log_mbh_vs_KMAG_sph()
	#sani_vs_me()
	#BTratios_GS13_vs_me()
	#mag_tot_GS13_vs_me()
	#mag_tot_vs_sigma()
	
main()
