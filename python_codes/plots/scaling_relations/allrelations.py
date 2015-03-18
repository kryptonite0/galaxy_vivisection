import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy
from instruments import b_n
from instruments.linear_regression import colorline


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 12})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_auxiliary_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/auxiliary/'
path_scalrel_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/scaling_relations/'

earlytypes = ['E', 'E/S0', 'S0']
latetypes = ['S0/S', 'S']
mergers = ['merger']

def all_plots():
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.bar, \
		anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.sigma, \
		pysres.mag_sph_eq_moffat_comb, pysres.mag_tot_eq_moffat_comb, \
		pysres.log_n_maj_moffat_comb, \
		pysres.mu_e_maj_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		errV.perr_log_n, errV.merr_log_n, \
		errV.perr_mu_e, errV.merr_mu_e, \
		errV.perr_mu_0, errV.merr_mu_0 \
        	FROM Ancillary AS anc \
        	JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
        	WHERE anc.fit1D_done = 1;'

	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	simplemorphtype = data[1]
	core = data[2].astype(np.int)
	bar = data[3].astype(np.int)
	
	mbh = data[4].astype(np.float)
	perr_mbh = data[5].astype(np.float)
	merr_mbh = data[6].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	
	sigma = data[7].astype(np.float)
	log_sigma = np.log10(sigma)
	err_log_sigma = sigma*[0.0] + np.log10(1.05)
	
	mag_sph = data[8].astype(np.float)
	mag_tot = data[9].astype(np.float)
	perr_mag_sph = data[12].astype(np.float)
	merr_mag_sph = data[13].astype(np.float)
	
	log_n = data[10].astype(np.float)
	n = 10**log_n
	perr_log_n = data[14].astype(np.float)
	merr_log_n = data[15].astype(np.float)
	
	mu_e = data[11].astype(np.float)
	perr_mu_e = data[16].astype(np.float)
	merr_mu_e = data[17].astype(np.float)
	
	# compute mu_0
	b = mu_e * [0.0]
	for i in range(len(b)):
		b[i] = b_n.computeb_n(n[i])
	
	mu_0 = mu_e - 2.5*b/np.log(10)
	perr_mu_0 = data[18].astype(np.float)
	merr_mu_0 = data[19].astype(np.float)

	## axis limits
	minsigma, maxsigma = 70, 450
	minmbh, maxmbh = 10**5.5, 10**10.9
	minmagsph, maxmagsph = -19.01, -28.99
	minmagtot, maxmagtot = -21.01, -28.99
	minn, maxn = 0.4, 15
	minmu0, maxmu0 = -10.1, -34.9
	
	############### panels ############### 
	
	## mbh vs sigma
	
	ax21 = plt.subplot(5,5,21) 
	
	outliers_mbh_sigma = ['n4889', 'n1374', 'n3842exp']
	out_mbh_sigma = [0]*sigma
	early_mbh_sigma = [0]*sigma
	late_mbh_sigma = [0]*sigma
	merger_mbh_sigma = [0]*sigma
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_mbh_sigma:
			out_mbh_sigma[i] = 1
		if gal_id[i] not in outliers_mbh_sigma and simplemorphtype[i] in earlytypes and sigma[i]>0:
			early_mbh_sigma[i] = 1	
		if gal_id[i] not in outliers_mbh_sigma and simplemorphtype[i] in latetypes and sigma[i]>0:
			late_mbh_sigma[i] = 1	
		if gal_id[i] not in outliers_mbh_sigma and simplemorphtype[i] in mergers and sigma[i]>0:
			merger_mbh_sigma[i] = 1	
			
	print 'mbh vs sigma'
	print 'all'
	print 'Sample size =', len(log_sigma[out_mbh_sigma==0])
        A,B,Aerr,Berr,covAB=bces.bces(log_sigma[out_mbh_sigma==0]-np.average(log_sigma[out_mbh_sigma==0]),
		err_log_sigma[out_mbh_sigma==0],
        	log_mbh[out_mbh_sigma==0],0.5*(merr_log_mbh[out_mbh_sigma==0] + perr_log_mbh[out_mbh_sigma==0]),log_sigma[out_mbh_sigma==0]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_sigma[out_mbh_sigma==0])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,1000,10)
	yy = (A[2]*logxx + B[2])
        ax21.plot(10**(logxx+np.average(log_sigma[out_mbh_sigma==0])),10**yy, color='k', ls='-', linewidth=1.5)
		
	ax21.scatter(sigma[early_mbh_sigma==1], mbh[early_mbh_sigma==1], marker='o', color='r', s=10)
	ax21.scatter(sigma[late_mbh_sigma==1], mbh[late_mbh_sigma==1], marker='o', color='b', s=10)
	ax21.scatter(sigma[merger_mbh_sigma==1], mbh[merger_mbh_sigma==1], marker='d', color='k', s=10)
	ax21.scatter(sigma[out_mbh_sigma==1], mbh[out_mbh_sigma==1], marker='x', color='k', s=10)

	ax21.axis([minsigma,maxsigma,minmbh,maxmbh])
	ax21.set_xscale('log')
	ax21.set_yscale('log')
	ax21.set_xlabel(r'$\sigma$', labelpad=5)
	ax21.set_ylabel(r'$M_{\rm BH}$', labelpad=5)
	
	## mbh vs magsph
	
	ax22 = plt.subplot(5,5,22) 
	
	outliers_mbh_mag_sph = ['n4889', 'n1374', 'n3842exp']
	out_mbh_mag_sph = [0]*mag_sph
	early_mbh_mag_sph = [0]*mag_sph
	late_mbh_mag_sph = [0]*mag_sph
	merger_mbh_mag_sph = [0]*mag_sph
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_mbh_mag_sph:
			out_mbh_mag_sph[i] = 1
		if gal_id[i] not in outliers_mbh_mag_sph and simplemorphtype[i] in earlytypes:
			early_mbh_mag_sph[i] = 1	
		if gal_id[i] not in outliers_mbh_mag_sph and simplemorphtype[i] in latetypes:
			late_mbh_mag_sph[i] = 1	
		if gal_id[i] not in outliers_mbh_mag_sph and simplemorphtype[i] in mergers:
			merger_mbh_mag_sph[i] = 1	
			
	print 'mbh vs magsph'
	print 'early'
	print 'Sample size =', len(mag_sph[early_mbh_mag_sph==1])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[early_mbh_mag_sph==1]-np.average(mag_sph[early_mbh_mag_sph==1]),
		0.5*(perr_mag_sph[early_mbh_mag_sph==1]+merr_mag_sph[early_mbh_mag_sph==1]),
        	log_mbh[early_mbh_mag_sph==1],0.5*(merr_log_mbh[early_mbh_mag_sph==1] + perr_log_mbh[early_mbh_mag_sph==1]),log_sigma[early_mbh_mag_sph==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[early_mbh_mag_sph==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax22.plot((logxx+np.average(mag_sph[early_mbh_mag_sph==1])),10**yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(mag_sph[late_mbh_mag_sph==1])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[late_mbh_mag_sph==1]-np.average(mag_sph[late_mbh_mag_sph==1]),
		0.5*(perr_mag_sph[late_mbh_mag_sph==1]+merr_mag_sph[late_mbh_mag_sph==1]),
        	log_mbh[late_mbh_mag_sph==1],0.5*(merr_log_mbh[late_mbh_mag_sph==1] + perr_log_mbh[late_mbh_mag_sph==1]),mag_sph[late_mbh_mag_sph==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[late_mbh_mag_sph==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax22.plot((logxx+np.average(mag_sph[late_mbh_mag_sph==1])),10**yy, color='b', ls='-', linewidth=1.5)
	
	ax22.scatter(mag_sph[early_mbh_mag_sph==1], mbh[early_mbh_mag_sph==1], marker='o', color='r', s=10)
	ax22.scatter(mag_sph[late_mbh_mag_sph==1], mbh[late_mbh_mag_sph==1], marker='o', color='b', s=10)
	ax21.scatter(mag_sph[merger_mbh_sigma==1], mbh[merger_mbh_sigma==1], marker='d', color='k', s=10)
	ax22.scatter(mag_sph[out_mbh_mag_sph==1], mbh[out_mbh_mag_sph==1], marker='x', color='k', s=10)

	ax22.axis([minmagsph,maxmagsph,minmbh,maxmbh])
	ax22.set_yscale('log')
	ax22.set_yticklabels([])
	ax22.set_xlabel(r'$MAG_{\rm sph}$', labelpad=5)
	
	## mbh vs magtot
	
	ax23 = plt.subplot(5,5,23) 
	
	ax23.axis([minmagtot,maxmagtot,minmbh,maxmbh])
	ax23.scatter(mag_tot, mbh, marker='o', color='k', s=10)
	ax23.scatter(mag_tot[gal_id=='n4889'], mbh[gal_id=='n4889'], marker='o', color='r', s=10)
	ax23.scatter(mag_tot[gal_id=='n1374'], mbh[gal_id=='n1374'], marker='o', color='r', s=10)
	ax23.scatter(mag_tot[gal_id=='n3842exp'], mbh[gal_id=='n3842exp'], marker='o', color='r', s=10)
	ax23.set_yscale('log')
	ax23.set_yticklabels([])
	ax23.set_xlabel(r'$MAG_{\rm tot}$', labelpad=5)
	
	## mbh vs n
	
	ax24 = plt.subplot(5,5,24) 
	
	outliers_mbh_n = ['n4889', 'n1374', 'n3842exp', 'n3998', 'n0524', 'n3377', 'n4697']
	out_mbh_n = [0]*n
	early_mbh_n = [0]*n
	late_mbh_n = [0]*n
	merger_mbh_n = [0]*n
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_mbh_n:
			out_mbh_n[i] = 1
		if gal_id[i] not in outliers_mbh_n and simplemorphtype[i] in earlytypes:
			early_mbh_n[i] = 1	
		if gal_id[i] not in outliers_mbh_n and simplemorphtype[i] in latetypes:
			late_mbh_n[i] = 1
		if gal_id[i] not in outliers_mbh_n and simplemorphtype[i] in mergers:
			merger_mbh_n[i] = 1	
	
	print 'mbh vs n'
	print 'early'
	print 'Sample size =', len(log_n[early_mbh_n==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[early_mbh_n==1]-np.average(log_n[early_mbh_n==1]),
		0.5*(perr_log_n[early_mbh_n==1]+merr_log_n[early_mbh_n==1]),
        	log_mbh[early_mbh_n==1],0.5*(merr_log_mbh[early_mbh_n==1] + perr_log_mbh[early_mbh_n==1]),log_sigma[early_mbh_n==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[early_mbh_n==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax24.plot(10**(logxx+np.average(log_n[early_mbh_n==1])),10**yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(log_n[late_mbh_n==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[late_mbh_n==1]-np.average(log_n[late_mbh_n==1]),
		0.5*(perr_log_n[late_mbh_n==1]+merr_log_n[late_mbh_n==1]),
        	log_mbh[late_mbh_n==1],0.5*(merr_log_mbh[late_mbh_n==1] + perr_log_mbh[late_mbh_n==1]),log_n[late_mbh_n==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[late_mbh_n==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax24.plot(10**(logxx+np.average(log_n[late_mbh_n==1])),10**yy, color='b', ls='-', linewidth=1.5)
			
	ax24.scatter(n[early_mbh_n==1], mbh[early_mbh_n==1], marker='o', color='r', s=10)
	ax24.scatter(n[late_mbh_n==1], mbh[late_mbh_n==1], marker='o', color='b', s=10)
	ax24.scatter(n[merger_mbh_n==1], mbh[merger_mbh_n==1], marker='d', color='k', s=10)
	ax24.scatter(n[out_mbh_n==1], mbh[out_mbh_n==1], marker='x', color='k', s=10)
		
	ax24.axis([minn,maxn,minmbh,maxmbh])
	ax24.set_xscale('log')
	ax24.set_yscale('log')
	ax24.set_yticklabels([])
	ax24.set_xlabel(r'$n_{\rm sph}$', labelpad=5)
	
	## mbh vs mu_0
	
	ax25 = plt.subplot(5,5,25) 
	
	outliers_mbh_mu0 = ['n4889', 'n1374', 'n3842exp', 'm31']
	out_mbh_mu0 = [0]*mu_0
	early_mbh_mu0 = [0]*mu_0
	late_mbh_mu0 = [0]*mu_0
	merger_mbh_mu0 = [0]*mu_0
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_mbh_mu0:
			out_mbh_mu0[i] = 1
		if gal_id[i] not in outliers_mbh_mu0 and simplemorphtype[i] in earlytypes:
			early_mbh_mu0[i] = 1	
		if gal_id[i] not in outliers_mbh_mu0 and simplemorphtype[i] in latetypes:
			late_mbh_mu0[i] = 1
		if gal_id[i] not in outliers_mbh_mu0 and simplemorphtype[i] in mergers:
			merger_mbh_mu0[i] = 1	
	
	print 'mbh vs mu0'
	print 'early'
	print 'Sample size =', len(mu_0[early_mbh_mu0==1])
        A,B,Aerr,Berr,covAB=bces.bces(mu_0[early_mbh_mu0==1]-np.average(mu_0[early_mbh_mu0==1]),
		0.5*(perr_mu_0[early_mbh_mu0==1]+merr_mu_0[early_mbh_mu0==1]),
        	log_mbh[early_mbh_mu0==1],0.5*(merr_log_mbh[early_mbh_mu0==1] + perr_log_mbh[early_mbh_mu0==1]),mu_0[early_mbh_mu0==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mu_0[early_mbh_mu0==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax25.plot((logxx+np.average(mu_0[early_mbh_mu0==1])),10**yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(mu_0[late_mbh_mu0==1])
        A,B,Aerr,Berr,covAB=bces.bces(mu_0[late_mbh_mu0==1]-np.average(mu_0[late_mbh_mu0==1]),
		0.5*(perr_mu_0[late_mbh_mu0==1]+merr_mu_0[late_mbh_mu0==1]),
        	log_mbh[late_mbh_mu0==1],0.5*(merr_log_mbh[late_mbh_mu0==1] + perr_log_mbh[late_mbh_mu0==1]),mu_0[late_mbh_mu0==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mu_0[late_mbh_mu0==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax25.plot((logxx+np.average(mu_0[late_mbh_mu0==1])),10**yy, color='b', ls='-', linewidth=1.5)
	
	ax25.scatter(mu_0[early_mbh_mu0==1], mbh[early_mbh_mu0==1], marker='o', color='r', s=10)
	ax25.scatter(mu_0[late_mbh_mu0==1], mbh[late_mbh_mu0==1], marker='o', color='b', s=10)
	ax25.scatter(mu_0[merger_mbh_n==1], mbh[merger_mbh_n==1], marker='d', color='k', s=10)
	ax25.scatter(mu_0[out_mbh_mu0==1], mbh[out_mbh_mu0==1], marker='x', color='k', s=10)
	
	ax25.axis([minmu0,maxmu0,minmbh,maxmbh])
	ax25.set_yscale('log')
	ax25.set_yticklabels([])
	ax25.set_xlabel(r'$\mu_{\rm 0,sph}$', labelpad=5)
	
	## sigma vs magsph
	
	ax17 = plt.subplot(5,5,17) 
	
	early_sigma_magsph = [0]*sigma
	late_sigma_magsph = [0]*sigma
	merger_sigma_magsph = [0]*sigma
	for i in range(len(gal_id)):
		if simplemorphtype[i] in earlytypes and sigma[i]>0:
			early_sigma_magsph[i] = 1	
		if simplemorphtype[i] in latetypes and sigma[i]>0:
			late_sigma_magsph[i] = 1
		if simplemorphtype[i] in mergers and sigma[i]>0:
			merger_sigma_magsph[i] = 1	
	
	print 'sigma vs magsph'
	print 'early'
	print 'Sample size =', len(mag_sph[early_sigma_magsph==1])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[early_sigma_magsph==1]-np.average(mag_sph[early_sigma_magsph==1]),
		0.5*(perr_mag_sph[early_sigma_magsph==1]+merr_mag_sph[early_sigma_magsph==1]),
        	log_sigma[early_sigma_magsph==1],err_log_sigma[early_sigma_magsph==1],mag_sph[early_sigma_magsph==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[early_sigma_magsph==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax17.plot((logxx+np.average(mag_sph[early_sigma_magsph==1])),10**yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(mag_sph[late_sigma_magsph==1])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[late_sigma_magsph==1]-np.average(mag_sph[late_sigma_magsph==1]),
		0.5*(perr_mag_sph[late_sigma_magsph==1]+merr_mag_sph[late_sigma_magsph==1]),
        	log_sigma[late_sigma_magsph==1],err_log_sigma[late_sigma_magsph==1],log_sigma[late_sigma_magsph==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[late_sigma_magsph==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
        yy = (A[2]*logxx + B[2])
        ax17.plot((logxx+np.average(mag_sph[late_sigma_magsph==1])),10**yy, color='b', ls='-', linewidth=1.5)

	ax17.scatter(mag_sph[early_sigma_magsph==1], sigma[early_sigma_magsph==1], marker='o', color='r', s=10)
	ax17.scatter(mag_sph[late_sigma_magsph==1], sigma[late_sigma_magsph==1], marker='o', color='b', s=10)
	ax17.scatter(mag_sph[merger_sigma_magsph==1], sigma[merger_sigma_magsph==1], marker='d', color='k', s=10)

	ax17.axis([minmagsph,maxmagsph,minsigma,maxsigma])
	ax17.set_yscale('log')
	ax17.set_xticklabels([])
	ax17.set_ylabel(r'$\sigma$', labelpad=5)
	
	## sigma vs magtot
	
	ax18 = plt.subplot(5,5,18) 
	
	ax18.axis([minmagtot,maxmagtot,minsigma,maxsigma])
	ax18.scatter(mag_tot, sigma, marker='o', color='k', s=10)
	ax18.set_yscale('log')
	ax18.set_xticklabels([])
	ax18.set_yticklabels([])
	
	## sigma vs n
	
	ax19 = plt.subplot(5,5,19) 
	
	outliers_sigma_n = ['n3998', 'n0524', 'n3377', 'n4697']
	out_sigma_n = [0]*n
	early_sigma_n = [0]*n
	late_sigma_n = [0]*n
	merger_sigma_n = [0]*n
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_sigma_n:
			out_sigma_n[i] = 1
		if gal_id[i] not in outliers_sigma_n and simplemorphtype[i] in earlytypes:
			early_sigma_n[i] = 1	
		if gal_id[i] not in outliers_sigma_n and simplemorphtype[i] in latetypes:
			late_sigma_n[i] = 1
		if gal_id[i] not in outliers_sigma_n and simplemorphtype[i] in mergers:
			merger_sigma_n[i] = 1	
	
	ax19.scatter(n, sigma, marker='o', color='k', s=10)
	ax19.scatter(n[early_sigma_n==1], sigma[early_sigma_n==1], marker='o', color='r', s=10)
	ax19.scatter(n[late_sigma_n==1], sigma[late_sigma_n==1], marker='o', color='b', s=10)
	ax19.scatter(n[merger_sigma_n==1], sigma[merger_sigma_n==1], marker='d', color='k', s=10)
	
	ax19.axis([minn,maxn,minsigma,maxsigma])
	ax19.set_xscale('log')
	ax19.set_yscale('log')
	ax19.set_xticklabels([])
	ax19.set_yticklabels([])
	
	## sigma vs mu_0
	
	ax20 = plt.subplot(5,5,20) 
	
	outliers_sigma_mu0 = ['m31']
	out_sigma_mu0 = [0]*mu_0
	early_sigma_mu0 = [0]*mu_0
	late_sigma_mu0 = [0]*mu_0
	merger_sigma_mu0 = [0]*mu_0
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_sigma_mu0:
			out_sigma_mu0[i] = 1
		if gal_id[i] not in outliers_sigma_mu0 and simplemorphtype[i] in earlytypes and sigma[i]>0:
			early_sigma_mu0[i] = 1	
		if gal_id[i] not in outliers_sigma_mu0 and simplemorphtype[i] in latetypes and sigma[i]>0:
			late_sigma_mu0[i] = 1
		if gal_id[i] not in outliers_sigma_mu0 and simplemorphtype[i] in mergers and sigma[i]>0:
			merger_sigma_mu0[i] = 1	
	
	print 'sigma vs mu0'
	print 'early'
	print 'Sample size =', len(mu_0[early_sigma_mu0==1])
        A,B,Aerr,Berr,covAB=bces.bces(mu_0[early_sigma_mu0==1]-np.average(mu_0[early_sigma_mu0==1]),
		0.5*(perr_mu_0[early_sigma_mu0==1]+merr_mu_0[early_sigma_mu0==1]),
        	log_sigma[early_sigma_mu0==1],err_log_sigma[early_sigma_mu0==1],mu_0[early_sigma_mu0==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mu_0[early_sigma_mu0==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax20.plot((logxx+np.average(mu_0[early_sigma_mu0==1])),10**yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(mu_0[late_sigma_mu0==1])
        A,B,Aerr,Berr,covAB=bces.bces(mu_0[late_sigma_mu0==1]-np.average(mu_0[late_sigma_mu0==1]),
		0.5*(perr_mu_0[late_sigma_mu0==1]+merr_mu_0[late_sigma_mu0==1]),
        	log_sigma[late_sigma_mu0==1],err_log_sigma[late_sigma_mu0==1],mu_0[late_sigma_mu0==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mu_0[late_sigma_mu0==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax20.plot((logxx+np.average(mu_0[late_sigma_mu0==1])),10**yy, color='b', ls='-', linewidth=1.5)
	
	ax20.scatter(mu_0[early_sigma_mu0==1], sigma[early_sigma_mu0==1], marker='o', color='r', s=10)
	ax20.scatter(mu_0[late_sigma_mu0==1], sigma[late_sigma_mu0==1], marker='o', color='b', s=10)
	ax20.scatter(mu_0[merger_sigma_n==1], sigma[merger_sigma_n==1], marker='d', color='k', s=10)
	ax20.scatter(mu_0[out_sigma_mu0==1], sigma[out_sigma_mu0==1], marker='x', color='k', s=10)
	
	ax20.axis([minmu0,maxmu0,minsigma,maxsigma])
	ax20.set_yscale('log')
	ax20.set_xticklabels([])
	ax20.set_yticklabels([])
	
	## magsph vs magtot
	
	ax13 = plt.subplot(5,5,13) 
	
	ax13.axis([minmagtot,maxmagtot,minmagsph,maxmagsph])
	ax13.scatter(mag_tot, mag_sph, marker='o', color='k', s=10)
	ax13.set_xticklabels([])
	ax13.set_ylabel(r'$MAG_{\rm sph}$', labelpad=5)
	
	## magsph vs n
	
	ax14 = plt.subplot(5,5,14) 
	
	outliers_magsph_n = ['n3998', 'n0524', 'n3377', 'n4697']
	out_magsph_n = [0]*n
	early_magsph_n = [0]*n
	late_magsph_n = [0]*n
	merger_magsph_n = [0]*n
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_magsph_n:
			out_magsph_n[i] = 1
		if gal_id[i] not in outliers_magsph_n and simplemorphtype[i] in earlytypes:
			early_magsph_n[i] = 1	
		if gal_id[i] not in outliers_magsph_n and simplemorphtype[i] in latetypes:
			late_magsph_n[i] = 1
		if gal_id[i] not in outliers_magsph_n and simplemorphtype[i] in mergers:
			merger_magsph_n[i] = 1	
	
	print 'magsph vs n'
	print 'early'
	print 'Sample size =', len(log_n[early_magsph_n==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[early_magsph_n==1]-np.average(log_n[early_magsph_n==1]),
		0.5*(perr_log_n[early_magsph_n==1]+merr_log_n[early_magsph_n==1]),
        	mag_sph[early_magsph_n==1],0.5*(merr_mag_sph[early_magsph_n==1] + perr_mag_sph[early_magsph_n==1]),n[early_magsph_n==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[early_magsph_n==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax14.plot(10**(logxx+np.average(log_n[early_magsph_n==1])),yy, color='r', ls='-', linewidth=1.5)
		
	print 'late'
	print 'Sample size =', len(log_n[late_magsph_n==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[late_magsph_n==1]-np.average(log_n[late_magsph_n==1]),
		0.5*(perr_log_n[late_magsph_n==1]+merr_log_n[late_magsph_n==1]),
        	mag_sph[late_magsph_n==1],0.5*(merr_mag_sph[late_magsph_n==1] + perr_mag_sph[late_magsph_n==1]),n[late_magsph_n==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[late_magsph_n==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-100,100,10)
	yy = (A[2]*logxx + B[2])
        ax14.plot(10**(logxx+np.average(log_n[late_magsph_n==1])),yy, color='b', ls='-', linewidth=1.5)
	
	ax14.scatter(n[early_magsph_n==1], mag_sph[early_magsph_n==1], marker='o', color='r', s=10)
	ax14.scatter(n[late_magsph_n==1], mag_sph[late_magsph_n==1], marker='o', color='b', s=10)
	ax14.scatter(n[merger_magsph_n==1], mag_sph[merger_magsph_n==1], marker='d', color='k', s=10)
	ax14.scatter(n[out_magsph_n==1], mag_sph[out_magsph_n==1], marker='x', color='k', s=10)
	
	ax14.axis([minn,maxn,minmagsph,maxmagsph])
	ax14.set_xscale('log')
	ax14.set_xticklabels([])
	ax14.set_yticklabels([])
	
	## magsph vs mu0
	
	ax15 = plt.subplot(5,5,15) 
	
	outliers_magsph_mu0 = ['n3998', 'n0524', 'n3377', 'n4697']
	out_magsph_mu0 = [0]*mu_0
	early_magsph_mu0 = [0]*mu_0
	late_magsph_mu0 = [0]*mu_0
	merger_magsph_mu0 = [0]*mu_0
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_magsph_mu0:
			out_magsph_mu0[i] = 1
		if gal_id[i] not in outliers_magsph_mu0 and simplemorphtype[i] in earlytypes:
			early_magsph_mu0[i] = 1	
		if gal_id[i] not in outliers_magsph_mu0 and simplemorphtype[i] in latetypes:
			late_magsph_mu0[i] = 1
		if gal_id[i] not in outliers_magsph_mu0 and simplemorphtype[i] in mergers:
			merger_magsph_mu0[i] = 1	
		
	ax15.scatter(mu_0[early_magsph_mu0==1], mag_sph[early_magsph_mu0==1], marker='o', color='r', s=10)
	ax15.scatter(mu_0[late_magsph_mu0==1], mag_sph[late_magsph_mu0==1], marker='o', color='b', s=10)
	ax15.scatter(mu_0[merger_magsph_mu0==1], mag_sph[merger_magsph_mu0==1], marker='d', color='k', s=10)
	ax15.scatter(mu_0[out_magsph_mu0==1], mag_sph[out_magsph_mu0==1], marker='x', color='k', s=10)
	
	ax15.axis([minmu0,maxmu0,minmagsph,maxmagsph])
	ax15.set_xticklabels([])
	ax15.set_yticklabels([])
	
	## magtot vs n
	
	ax9 = plt.subplot(5,5,9) 
	
	ax9.axis([minn,maxn,minmagtot,maxmagtot])
	ax9.scatter(n, mag_tot, marker='o', color='k', s=10)
	ax9.scatter(n[gal_id=='n3998'], mag_tot[gal_id=='n3998'], marker='o', color='g', s=10)
	ax9.scatter(n[gal_id=='n0524'], mag_tot[gal_id=='n0524'], marker='o', color='g', s=10)
	ax9.scatter(n[gal_id=='n3377'], mag_tot[gal_id=='n3377'], marker='o', color='g', s=10)
	ax9.scatter(n[gal_id=='n4697'], mag_tot[gal_id=='n4697'], marker='o', color='g', s=10)
	ax9.set_xscale('log')
	ax9.set_xticklabels([])
	ax9.set_ylabel(r'$MAG_{\rm tot}$', labelpad=5)
	
	## magtot vs mu0
	
	ax10 = plt.subplot(5,5,10) 
	
	ax10.axis([minmu0,maxmu0,minmagtot,maxmagtot])
	ax10.scatter(mu_0, mag_tot, marker='o', color='k', s=10)
	ax10.scatter(mu_0[gal_id=='m31'], mag_tot[gal_id=='m31'], marker='o', color='pink', s=10)
	ax10.set_xticklabels([])
	ax10.set_xticklabels([])
	ax10.set_yticklabels([])
	
	## n vs mu0
	
	ax5 = plt.subplot(5,5,5) 
	
	outliers_n_mu0 = ['m31']
	out_n_mu0 = [0]*mu_0
	early_n_mu0 = [0]*mu_0
	late_n_mu0 = [0]*mu_0
	merger_n_mu0 = [0]*mu_0
	for i in range(len(gal_id)):
		if gal_id[i] in outliers_n_mu0:
			out_n_mu0[i] = 1
		if gal_id[i] not in outliers_n_mu0 and simplemorphtype[i] in earlytypes:
			early_n_mu0[i] = 1	
		if gal_id[i] not in outliers_n_mu0 and simplemorphtype[i] in latetypes:
			late_n_mu0[i] = 1
		if gal_id[i] not in outliers_n_mu0 and simplemorphtype[i] in mergers:
			merger_n_mu0[i] = 1	
		
	ax5.scatter(mu_0[early_n_mu0==1], n[early_n_mu0==1], marker='o', color='r', s=10)
	ax5.scatter(mu_0[late_n_mu0==1], n[late_n_mu0==1], marker='o', color='b', s=10)
	ax5.scatter(mu_0[merger_n_mu0==1], n[merger_n_mu0==1], marker='d', color='k', s=10)
	ax5.scatter(mu_0[out_n_mu0==1], n[out_n_mu0==1], marker='x', color='k', s=10)
	
	ax5.axis([minmu0,maxmu0,minn,maxn])
	ax5.set_yscale('log')
	ax5.set_xticklabels([])
	ax5.set_ylabel(r'$n_{\rm sph}$', labelpad=5)
	

	plt.show()
	#plt.savefig(path_scalrel_plots + 'all_plots.pdf', format='pdf', dpi=1000)
	
all_plots()	
	
	
