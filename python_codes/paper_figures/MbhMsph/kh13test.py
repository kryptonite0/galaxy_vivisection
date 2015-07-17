import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.image as image

from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy
from instruments import b_n
from instruments.linear_regression import colorline
from instruments.linear_regression import akaike
from instruments.linear_regression import absolutescatter
#from instruments.linear_regression import lnr
from instruments.markers import markers

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 26})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/MbhMsph/images/'


def mbh_vs_mass_sph():
	
	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, kh13.mass_BH, kh13.perr_mass_BH, kh13.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype, \
		pysres.log_n_eq_moffat_comb, \
		col.color  \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		JOIN Colors as col ON anc.gal_id = col.gal_id \
		JOIN KormendyHo2013 as kh13 ON anc.gal_id = kh13.gal_id \
		WHERE anc.fit1D_done = 1 \
		AND anc.gal_id != "ic2560" \
		AND anc.gal_id != "n2778" ;'

	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
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
	simplemorphtype = data[9]
	log_n = data[10].astype(np.float)
	n = 10**log_n
	color = data[11].astype(np.float)
	
	#print mbh, merr_mbh
	print len(log_mbh)
	
	log_ML = 3.98*color+0.13 # meidt+2014
	ML = 10**log_ML
	
	mass_sph = ML*10**(-0.4*(mag_sph-3.25))
	merr_mass_sph = -ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) + mass_sph
	perr_mass_sph = +ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) - mass_sph
	
	log_mass_sph = np.log10(mass_sph)
	perr_log_mass_sph = np.log10(1+perr_mass_sph/mass_sph)
	merr_log_mass_sph = -np.log10(1-merr_mass_sph/mass_sph)
		
        ELL_core = ELL*[0]
        ELL_sersic = ELL*[0]
        BUL_core = ELL*[0]
        BUL_sersic = ELL*[0]
       
        earlytype = ELL*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = core*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
	
        for i in range(len(ELL)):
        	if simplemorphtype[i]=='E' and core[i]==1 and gal_id[i] not in outliers:
        		ELL_core[i] = 1
        	elif simplemorphtype[i]=='E' and core[i]==0 and gal_id[i] not in outliers:
        		ELL_sersic[i] = 1
        	elif bulge[i]==1 and core[i]==1 and gal_id[i] not in outliers:
        		BUL_core[i] = 1
        	elif bulge[i]==1 and core[i]==0 and gal_id[i] not in outliers:
        		BUL_sersic[i] = 1
        
	morph_coreList = []
        for i in range(len(simplemorphtype)):
        	if gal_id[i] not in outliers:
			morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
		else:
			morph_coreList.append('outlier')	
        morph_core = np.asarray(morph_coreList)
	       
        ### produce .dat file for paper with Ewan (19 June 2015)
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/tables/mbh_vs_mass_sph.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# galaxy    logMassSph     +err_logMassSph     -err_logMassSph     average_err_logMassSph     logMassBH     +err_logMassBH    -err_logMassBH    average_err_logMassBH     core   type   \n')
        for a,b,c,d,e,f,g,h,i,j,k in zip(gal_id, log_mass_sph, perr_log_mass_sph, merr_log_mass_sph, 
		0.5*(perr_log_mass_sph + merr_log_mass_sph),
        	log_mbh, perr_log_mbh, merr_log_mbh, 0.5*(merr_log_mbh + perr_log_mbh), 
		core, simplemorphtype):
        	datfile.write(str(a) + '    ' + str(b) + '    ' + str(c) + '    ' + str(d) + '    ')
		datfile.write(str(e) + '    ' + str(f) + '    ' + str(g) + '    ' + str(h) + '    ') 
		datfile.write(str(i) + '    ' + str(j) + '    ' + str(k) + '\n')
        datfile.close()

        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

	for x0,y0 in zip(mass_sph[simplemorphtype=='E'], mbh[simplemorphtype=='E']):
		markers.elliptical(ax, 'red', np.log10(x0), np.log10(y0), 0.08)
	
 	for x0,y0 in zip(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
		markers.lenticular(ax, 'red', np.log10(x0), np.log10(y0), 0.08)
       
	for x0,y0 in zip(mass_sph[simplemorphtype=='S0'], mbh[simplemorphtype=='S0']):
		markers.lenticular(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.08)
	
        for x0,y0 in zip(mass_sph[simplemorphtype=='S0/Sp'], mbh[simplemorphtype=='S0/Sp']):
		markers.spiral(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.04)
			
        for x0,y0 in zip(mass_sph[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
		markers.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.04)
		
        for x0,y0 in zip(mass_sph[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
		markers.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.04)

	ax.scatter(mass_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
			
	
	print 'early'
	print 'n', len(log_mass_sph[earlytype==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),
		0.5*(perr_log_mass_sph[earlytype==1] + merr_log_mass_sph[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),log_mass_sph[earlytype==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph[earlytype==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '   B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '   B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1])
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '   B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy, color='r', ls='--', linewidth=2.)
	#colorline.colorline(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**yy, cmap=green_red)
	
	##### calculates 1sigma uncertainty band
	yy_1 = ((A[2]+Aerr[2])*(logxx) + (B[2]+Berr[2]))
	yy_2 = ((A[2]-Aerr[2])*(logxx) + (B[2]+Berr[2]))
	yy_3 = ((A[2]+Aerr[2])*(logxx) + (B[2]-Berr[2]))
	yy_4 = ((A[2]-Aerr[2])*(logxx) + (B[2]-Berr[2]))
	yy_up = yy_1*[0.0]
	for i in range(len(yy_1)):
		if yy_1[i] > yy_2[i]:
			yy_up[i] = yy_1[i]
		elif yy_1[i] <= yy_2[i]:
			yy_up[i] = yy_2[i]	
	yy_up = np.asarray(yy_up)
	yy_lo = yy_1*[0.0]
	for i in range(len(yy_3)):
		if yy_3[i] < yy_4[i]:
			yy_lo[i] = yy_3[i]
		elif yy_3[i] >= yy_4[i]:
			yy_lo[i] = yy_4[i]	
	yy_lo = np.asarray(yy_lo)
				
	ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.3, facecolor='r')
	
        print 'sersic bul of spi'
        print 'n', len(log_mass_sph[morph_core=='Sp_0'])
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[morph_core=='Sp_0']-np.average(log_mass_sph[morph_core=='Sp_0']),
		0.5*(perr_log_mass_sph[morph_core=='Sp_0'] + merr_log_mass_sph[morph_core=='Sp_0']),
        	log_mbh[morph_core=='Sp_0'],0.5*(merr_log_mbh[morph_core=='Sp_0'] + perr_log_mbh[morph_core=='Sp_0']),log_mass_sph[morph_core=='Sp_0']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph[morph_core=='Sp_0'])
        print
        #print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '   B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '   B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1])
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '   B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(10**(logxx+np.average(log_mass_sph[morph_core=='Sp_0'])),10**yy, color='b', ls='-', linewidth=2.)
	
	##### calculates 1sigma uncertainty band
	yy_1 = ((A[2]+Aerr[2])*(logxx) + (B[2]+Berr[2]))
	yy_2 = ((A[2]-Aerr[2])*(logxx) + (B[2]+Berr[2]))
	yy_3 = ((A[2]+Aerr[2])*(logxx) + (B[2]-Berr[2]))
	yy_4 = ((A[2]-Aerr[2])*(logxx) + (B[2]-Berr[2]))
	yy_up = yy_1*[0.0]
	for i in range(len(yy_1)):
		if yy_1[i] > yy_2[i]:
			yy_up[i] = yy_1[i]
		elif yy_1[i] <= yy_2[i]:
			yy_up[i] = yy_2[i]	
	yy_up = np.asarray(yy_up)
	yy_lo = yy_1*[0.0]
	for i in range(len(yy_3)):
		if yy_3[i] < yy_4[i]:
			yy_lo[i] = yy_3[i]
		elif yy_3[i] >= yy_4[i]:
			yy_lo[i] = yy_4[i]	
	yy_lo = np.asarray(yy_lo)
				
	ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='Sp_0'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='b')
	
	# legend
	markers.elliptical(ax, 'red', 8.85, 10.9, 0.08)
	ax.text(10**9.1, 10**10.75, 'E')
	markers.lenticular(ax, 'red', 8.85, 10.45, 0.08)
	ax.text(10**9.1, 10**10.3, 'E/S0')
	markers.lenticular(ax, 'darkorange', 8.85, 10., 0.08)
	ax.text(10**9.1, 10**9.85, 'S0')
	markers.spiral(ax, 'darkorange', 9.8, 10.9, 0.04)
	ax.text(10**10.05, 10**10.75, 'S0/Sp')
	markers.spiral(ax, 'blue', 9.8, 10.45, 0.04)
	ax.text(10**10.05, 10**10.3, 'Sp')
	ax.scatter([10**9.8], [10**10.], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
	ax.text(10**10.05, 10**9.85, 'merger')
	
 	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**8.6,10**12.3,10**5.3,10**11.2])
        plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mass_sph.pdf', format='pdf', dpi=1000)



def mbh_vs_mag_sph_psb():

	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, kh13.mass_BH, kh13.perr_mass_BH, kh13.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype, \
		pysres.log_n_eq_moffat_comb, \
		anc.pseudobulge_KH13 \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		JOIN KormendyHo2013 as kh13 ON anc.gal_id = kh13.gal_id \
		WHERE anc.fit1D_done = 1 \
		AND anc.gal_id != "ic2560" \
		AND anc.gal_id != "n2778" ;'
	
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
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
	simplemorphtype = data[9]
	log_n = data[10].astype(np.float)
	n = 10**log_n
	pseudobulge_KH13 = data[11].astype(np.int)
		
        ELL_core = ELL*[0]
        ELL_sersic = ELL*[0]
        BUL_core = ELL*[0]
        BUL_sersic = ELL*[0]
       
        earlytype = ELL*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = core*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
	
        for i in range(len(ELL)):
        	if simplemorphtype[i]=='E' and core[i]==1 and gal_id[i] not in outliers:
        		ELL_core[i] = 1
        	elif simplemorphtype[i]=='E' and core[i]==0 and gal_id[i] not in outliers:
        		ELL_sersic[i] = 1
        	elif bulge[i]==1 and core[i]==1 and gal_id[i] not in outliers:
        		BUL_core[i] = 1
        	elif bulge[i]==1 and core[i]==0 and gal_id[i] not in outliers:
        		BUL_sersic[i] = 1
        
	morph_coreList = []
        for i in range(len(simplemorphtype)):
        	if gal_id[i] not in outliers:
			morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
		else:
			morph_coreList.append('outlier')	
        morph_core = np.asarray(morph_coreList)
	
        
        fig, ax = plt.subplots()
	
	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

	ax.errorbar(mag_sph, mbh, xerr=[merr_mag_sph,perr_mag_sph], yerr=[merr_mbh,perr_mbh], 
		ecolor='gray', marker=' ', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False, **error_kwargs)
	
	cp = ax.scatter(mag_sph, mbh, marker='o', c=n, edgecolors='gray', s=200, **scatter_kwargs)
	ax.scatter(mag_sph[n<2], mbh[n<2], marker='s', facecolors='none', edgecolors='k', s=250)
	
	#ax.scatter([10**8.8], [10**10.7], marker='s', facecolors='none', edgecolors='k', s=300, **scatter_kwargs)
	#ax.scatter([10**8.8], [10**10.2], marker='*', facecolors='white', edgecolors='k', s=150, **scatter_kwargs)
	#ax.text(10**9.05,10**10.6, 'pseudobulge KH13')
	#ax.text(10**9.05,10**10.1, r'$n_{\rm sph}<2$')
	
	cbar = fig.colorbar(cp)
	#cbar.ax.set_ylabel(r'$n_{\rm sph}$', rotation=90)
	cbar.ax.set_ylabel(r'spheroid central radial concentration', rotation=90)
	
       #print 'BCES all'
       #print 'n', len(mag_sph)
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph-np.average(mag_sph),
       #	0.5*(perr_mag_sph + merr_mag_sph),
       #	log_mbh,0.5*(merr_log_mbh + perr_log_mbh),mag_sph*[0.0])
       ##absscat_0 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[0], A[0])
       ##absscat_1 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[2], A[2])
       ##absscat_3 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph)
       #print
       ##print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '   B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
       ##print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '   B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
       #print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '   B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot((logxx+np.average(mag_sph)),10**yy, color='k', ls='-', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY all'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph-np.average(mag_sph),
       #	0.5*(perr_mag_sph + merr_mag_sph),
       #	log_mbh,0.5*(merr_log_mbh + perr_log_mbh))
       #print '-----------------------------'
       #print '-----------------------------'
	
	###########################

	print 'BCES n>2'
	print 'n', len(mag_sph[n>2])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[n>2]-np.average(mag_sph[n>2]),
		0.5*(perr_mag_sph[n>2] + merr_mag_sph[n>2]),
        	log_mbh[n>2],0.5*(merr_log_mbh[n>2] + perr_log_mbh[n>2]),mag_sph[n>2]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(mag_sph[n>2]-np.average(mag_sph[n>2]), log_mbh[n>2], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(mag_sph[n>2]-np.average(mag_sph[n>2]), log_mbh[n>2], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(mag_sph[n>2]-np.average(mag_sph[n>2]), log_mbh[n>2], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(mag_sph[n>2]-np.average(mag_sph[n>2]), log_mbh[n>2], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[n>2])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '   B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
	print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '   B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '   B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
        
	slope_n_grt_2 = A[2]
	intercept_n_grt_2 = B[2]
	
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot((logxx+np.average(mag_sph[n>2])),10**yy, color='k', ls='-', linewidth=2.)
	
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY n>2'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[n>2]-np.average(mag_sph[n>2]),
       #	0.5*(perr_mag_sph[n>2] + merr_mag_sph[n>2]),
       #	log_mbh[n>2],0.5*(merr_log_mbh[n>2] + perr_log_mbh[n>2]))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #### produce .dat file
       #datfileName = '/Users/gsavorgnan/galaxy_vivisection/data/marconi_fit/mbh_vs_mag_sph_classbul.dat'
       #datfile = open(datfileName, 'w')
       #datfile.write('# MAGSph err_MAGSph logMassBH err_logMassBH \n')
       #for MAGSph, err_MAGSph, logMassBH, err_logMassBH in zip(mag_sph[n>2]-np.average(mag_sph[n>2]),
       #	0.5*(perr_mag_sph[n>2] + merr_mag_sph[n>2]),
       #	log_mbh[n>2],0.5*(merr_log_mbh[n>2] + perr_log_mbh[n>2])):
       #	datfile.write(str(MAGSph) + ' ' + str(err_MAGSph) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
       #datfile.close()

	###########################

       ## make inset
       #ins = plt.axes([.62, .25, .18, .2])
       #ins.axis([0.01,11,-1.99,1.99])
       #ins.xaxis.set_ticks([1,5,10])
       #ins.yaxis.set_ticks([-1,0,1])
       #ins.set_xlabel(r'$n_{\rm sph}$')
       #ins.set_ylabel(r'offset')
       #ins.plot([-1,15], [0,0], ls='-', c='k')
       #ins.plot([2,2], [-3,3], ls='--', c='k')
       #offset = log_mbh - (slope_n_grt_2*(mag_sph-np.average(mag_sph[n>2])) + intercept_n_grt_2)
       #ins.scatter(n, offset, marker='o', facecolor='k', edgecolor='k', s=15, **scatter_kwargs)
       ##ins.scatter(n[pseudobulge_KH13==1], offset[pseudobulge_KH13==1], marker='o', facecolor='b', edgecolor='gray', s=20, **scatter_kwargs)
	
	ax.set_yscale('log')
        ax.axis([-19.01,-27.99,10**5.1,10**11.2])
        xticks_labels = ['$-28$','','$-26$','','$-24$','','$-22$','','$-20$']
        ax.set_xticklabels(xticks_labels)
        #ax.set_xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=12)
        #ax.set_ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=12)
        ax.set_xlabel(r'spheroid absolute magnitude $\rm [mag]$', labelpad=12)
        ax.set_ylabel(r'black hole mass $\rm [M_\odot]$', labelpad=12)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mag_sph_psb.pdf', format='pdf', dpi=1000)
	#plt.savefig('mbh_vs_mag_sph_psb.pdf', format='pdf', dpi=1000)




def main():
	#mbh_vs_mass_sph()
	mbh_vs_mag_sph_psb()
	
	
main()	
