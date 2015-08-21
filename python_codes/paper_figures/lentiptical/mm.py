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

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/lentiptical/images/'

def mbh_vs_mass_sph():
	
	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype, \
		pysres.log_n_eq_moffat_comb, \
		col.color  \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		JOIN Colors as col ON anc.gal_id = col.gal_id \
		WHERE anc.fit1D_done = 1;'
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
		markers.elliptical(ax, 'k', np.log10(x0), np.log10(y0), 0.05)
	
 	#for x0,y0 in zip(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
	#	markers.lenticular(ax, 'red', np.log10(x0), np.log10(y0), 0.08)
        ax.scatter(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	
	mass_sph_n1277 = 2.69*10**11
	mass_sph_n1277_b = mass_sph_n1277/11.65*6	
	mass_sph_n1277_old = 2.88*10**10
	mbh_n1277 = 1.7*10**10
	ax.scatter([mass_sph_n1277], [mbh_n1277], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	#ax.scatter([mass_sph_n1277_b], [mbh_n1277], marker=r'$\star$', s=100, color='red', **scatter_kwargs)
	ax.scatter([mass_sph_n1277_old], [mbh_n1277], marker=r'$\star$', s=500, color='gray', **scatter_kwargs)
	ax.plot([mass_sph_n1277_old,mass_sph_n1277], [mbh_n1277,mbh_n1277], color='gray', lw=2, ls='--')
	
	mass_sph_n1271 = 9.17*10**10
	mass_sph_n1271_old = 5.4*10**10
	mbh_n1271 = 3*10**9
	ax.scatter([mass_sph_n1271], [mbh_n1271], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	#ax.plot([mass_sph_n1271_old,mass_sph_n1271], [mbh_n1271,mbh_n1271], color='blue', lw=3)
	
	mass_sph_m1216 = 2.27*10**11
	mbh_m1216 = 10**10 ## upper limit
	#ax.scatter([mass_sph_m1216], [mbh_m1216], marker='.', s=500, color='red', **scatter_kwargs)
	ax.scatter([mass_sph_m1216], [mbh_m1216/1.5], marker=r'$\downarrow$', s=500, color='red', **scatter_kwargs)
	
	mass_sph_n1332 = mass_sph[gal_id=='n1332']
	mass_sph_n1332_old = mass_sph_n1332/0.95*0.43
	mbh_n1332 = mbh[gal_id=='n1332']
	#ax.plot([mass_sph_n1332_old,mass_sph_n1332], [mbh_n1332,mbh_n1332], color='blue', lw=3)
	
	for x0,y0 in zip(mass_sph[simplemorphtype=='S0'], mbh[simplemorphtype=='S0']):
		markers.lenticular(ax, 'k', np.log10(x0), np.log10(y0), 0.05)
	
	x0 = mass_sph[gal_id=='n3115']
	y0 = mbh[gal_id=='n3115']
	ax.text(x0/1.5, 1.55*y0, 'N3115', size=12, color='red')
	
	x0 = mass_sph[gal_id=='n1332']
	y0 = mbh[gal_id=='n1332']
	ax.text(x0*1.3, y0/1.2, 'N1332', size=12, color='red')
	
	x0 = mass_sph_n1277
	y0 = mbh_n1277
	ax.text(x0/1.3, 1.6*y0, 'N1277', size=12, color='red')
	
	x0 = mass_sph_n1271
	y0 = mbh_n1271
	ax.text(x0/1.5, 1.6*y0, 'N1271', size=12, color='red')
	
	x0 = mass_sph_m1216
	y0 = mbh_m1216
	ax.text(x0/1.9, 0.7*y0, 'M1216', size=12, color='red')
	
	x0 = float(mass_sph[gal_id=='n4291'])
	y0 = float(mbh[gal_id=='n4291'])
	markers.elliptical(ax, 'green', np.log10(x0), np.log10(y0), 0.05)
	ax.text(x0/2, y0/1.6, 'N4291', size=12, color='green')
	
	x0 = float(mass_sph[gal_id=='n3998'])
	y0 = float(mbh[gal_id=='n3998'])
	#ax.text(x0/0.95, 1.3*y0, 'N3998', size=12, color='k')
	
        #for x0,y0 in zip(mass_sph[simplemorphtype=='S0/Sp'], mbh[simplemorphtype=='S0/Sp']):
	#	markers.spiral(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.04)
			
        #for x0,y0 in zip(mass_sph[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
	#	markers.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.04)
		
        #for x0,y0 in zip(mass_sph[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
	#	markers.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.04)

	#ax.scatter(mass_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
			
	
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
        ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy, color='k', ls='-', linewidth=2.)

	#epsilon = 0.43 # intrinsic scatter
	epsilon = 0.51 # tot rms scatter
	
	yy_1epsilon_up = (A[2]*(logxx) + B[2]) + epsilon
	yy_1epsilon_lo = (A[2]*(logxx) + B[2]) - epsilon
	ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_1epsilon_up, color='k', ls='--', linewidth=2.)
	ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_1epsilon_lo, color='k', ls='--', linewidth=2.)
	ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_1epsilon_lo,10**yy_1epsilon_up,alpha=0.3, facecolor='gray') 
		
	yy_3epsilon_up = (A[2]*(logxx) + B[2]) + 3*epsilon
	yy_3epsilon_lo = (A[2]*(logxx) + B[2]) - 3*epsilon
	ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_3epsilon_up, color='k', ls='--', linewidth=2.)
	ax.plot(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_3epsilon_lo, color='k', ls='--', linewidth=2.)
	ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])),10**yy_3epsilon_lo,10**yy_3epsilon_up,alpha=0.2, facecolor='gray') 
	
	ax.text(10**8.7, 10**6.87, r'$\pm 1 \sigma$', rotation=30, color='k', size=20)
	ax.text(10**8.7, 10**7.9, r'$\pm 3 \sigma$', rotation=30, color='k', size=20)
		
       ###### calculates 1sigma uncertainty band
       #yy_1 = ((A[2]+Aerr[2])*(logxx) + (B[2]+Berr[2]))
       #yy_2 = ((A[2]-Aerr[2])*(logxx) + (B[2]+Berr[2]))
       #yy_3 = ((A[2]+Aerr[2])*(logxx) + (B[2]-Berr[2]))
       #yy_4 = ((A[2]-Aerr[2])*(logxx) + (B[2]-Berr[2]))
       #yy_up = yy_1*[0.0]
       #for i in range(len(yy_1)):
       #	if yy_1[i] > yy_2[i]:
       #		yy_up[i] = yy_1[i]
       #	elif yy_1[i] <= yy_2[i]:
       #		yy_up[i] = yy_2[i]	
       #yy_up = np.asarray(yy_up)
       #yy_lo = yy_1*[0.0]
       #for i in range(len(yy_3)):
       #	if yy_3[i] < yy_4[i]:
       #		yy_lo[i] = yy_3[i]
       #	elif yy_3[i] >= yy_4[i]:
       #		yy_lo[i] = yy_4[i]	
       #yy_lo = np.asarray(yy_lo)
       #			
       #ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.3, facecolor='gray')
	
	# legend
	markers.elliptical(ax, 'k', 8.85, 10.9, 0.05)
	ax.text(10**9.1, 10**10.75, 'Elliptical')
	ax.scatter([10**8.85], [10**10.45], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	#markers.lenticular(ax, 'red', 8.85, 10.45, 0.08)
	ax.text(10**9.1, 10**10.3, 'Ellicular')
	markers.lenticular(ax, 'k', 8.85, 10., 0.05)
	ax.text(10**9.1, 10**9.85, 'Lenticular')
	#markers.spiral(ax, 'darkorange', 9.8, 10.9, 0.04)
	#ax.text(10**10.05, 10**10.75, 'S0/Sp')
	#markers.spiral(ax, 'blue', 9.8, 10.45, 0.04)
	#ax.text(10**10.05, 10**10.3, 'Sp')
	#ax.scatter([10**9.8], [10**10.], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
	#ax.text(10**10.05, 10**9.85, 'merger')
	
 	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**8.6,10**12.3,10**6.3,10**11.2])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
	plt.xlabel(r'Spheroid stellar mass $[M_\odot]$', labelpad=13)
	plt.ylabel(r'Black hole mass $[M_\odot]$', labelpad=13)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        #plt.show()
	plt.savefig(path_paper_figures + 'mm.pdf', format='pdf', dpi=1000)



mbh_vs_mass_sph()
