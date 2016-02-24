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
from instruments.markers import markers_n
from instruments.markers import markers_linlog

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 22})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/MbhN/images/'


def mag_sph_vs_n_maj():
	outliers = [u'n0524', u'n3998']
	#outliers = []
	
	drawfit = {'all' : True, 'E' : False, 'S0' : False, 'Sp' : True, 'bulge' : False, 'early' : True}
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, \
		physres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		physres.log_n_maj_moffat_comb, \
		errV.perr_log_n, errV.merr_log_n \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	simplemorphtype = data[1]
	core = data[2].astype(np.int)
	mag_sph = data[3].astype(np.float)
	perr_mag_sph = data[4].astype(np.float)
	merr_mag_sph = data[5].astype(np.float)
	log_n = data[6].astype(np.float)
	n = 10**log_n
	perr_log_n = data[7].astype(np.float)
	merr_log_n = data[8].astype(np.float)
	perr_n = (10**perr_log_n - 1)*n
	merr_n = (10**(-perr_log_n) - 1)*(-n)
	
	simplemorphtype[gal_id=='n0524'] = 'out'
	simplemorphtype[gal_id=='n3998'] = 'out'
	
	earlytype = mag_sph*[0]
        earlytype[simplemorphtype=='E'] = 1
        earlytype[simplemorphtype=='S0'] = 1
        earlytype[simplemorphtype=='E/S0'] = 1
	
	bulges = mag_sph*[0]
        bulges[simplemorphtype=='S0'] = 1
        bulges[simplemorphtype=='S0/Sp'] = 1
        bulges[simplemorphtype=='Sp'] = 1
	
	all = mag_sph*[0] + 1
	all[simplemorphtype=='out'] = 0
	#all[simplemorphtype=='merger'] = 0
	
	# build figure
	
	fig, (ax, ax_eb) = plt.subplots(2, sharex=True, sharey=True)
	
	scatter_kwargs = {"zorder":200}
	error_kwargs = {"lw":.5, "zorder":0}
	
	for x0,y0 in zip(n[simplemorphtype=='E'], mag_sph[simplemorphtype=='E']):
		markers_linlog.elliptical(ax, 'red', np.log10(x0), (y0), 0.035, 0.12)
	
        for x0,y0 in zip(n[simplemorphtype=='E/S0'], mag_sph[simplemorphtype=='E/S0']):
        	markers_linlog.lenticular(ax, 'red', np.log10(x0), (y0), 0.03, 0.12)
       
        for x0,y0 in zip(n[simplemorphtype=='S0'], mag_sph[simplemorphtype=='S0']):
        	markers_linlog.lenticular(ax, 'darkorange', np.log10(x0), (y0), 0.03, 0.12)
       
        for x0,y0 in zip(n[simplemorphtype=='S0/Sp'], mag_sph[simplemorphtype=='S0/Sp']):
        	markers_linlog.spiral(ax, 'darkorange', np.log10(x0), (y0), 0.02, 0.07)
        		
        for x0,y0 in zip(n[simplemorphtype=='Sp'], mag_sph[simplemorphtype=='Sp']):
        	markers_linlog.spiral(ax, 'blue', np.log10(x0), (y0), 0.02, 0.07)
        for x0,y0 in zip(n[simplemorphtype=='Sp'], mag_sph[simplemorphtype=='Sp']):
        	markers_linlog.spiral(ax, 'blue', np.log10(x0), (y0), 0.02, 0.07)
        	
        ax.scatter(n[simplemorphtype=='merger'], mag_sph[simplemorphtype=='merger'], marker=r'$\star$', s=200, color='k', **scatter_kwargs) 
	
	ax.scatter(n[gal_id=='n0524'], mag_sph[gal_id=='n0524'], marker='x', color='darkorange', s=100, linewidth=3)
	ax.scatter(n[gal_id=='n3998'], mag_sph[gal_id=='n3998'], marker='x', color='darkorange', s=100, linewidth=3)
	
	#expected_mag_sph = -23.94 -7.08*(log_n-0.51)
	#diff = mag_sph - expected_mag_sph
	#twosout = gal_id[diff>2*1.16]
	#threesout = gal_id[diff>3*1.16]
	#print twosout, threesout

	#######################################
	### linear regressions
	
        print 'all'
        print 'n', len(log_n[all==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[all==1]-np.average(log_n[all==1]),
        	0.5*(perr_log_n[all==1] + merr_log_n[all==1]),
        	mag_sph[all==1],0.5*(merr_mag_sph[all==1] + perr_mag_sph[all==1]),log_n[all==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), mag_sph[all==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), mag_sph[all==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), mag_sph[all==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), mag_sph[all==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[all==1])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['all'] == True:
		B_mfitexy = -23.90
		Berr_mfitexy = 0.13
		A_mfitexy = -7.04
		Aerr_mfitexy = 0.35
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[all==1])),yy, color='k', ls='dotted', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[all==1])),yy, color='k', ls='dotted', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[all==1])), yy_lo, yy_up, alpha=0.2, facecolor='k')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[all==1])), yy_lo, yy_up, alpha=0.2, facecolor='k')
              
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_all.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[all==1]-np.average(log_n[all==1]),
        	0.5*(perr_log_n[all==1] + merr_log_n[all==1]),
        	mag_sph[all==1],0.5*(merr_mag_sph[all==1] + perr_mag_sph[all==1])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
	
	
	#######################################
	
        print 'early'
        print 'n', len(log_n[earlytype==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
        	0.5*(perr_log_n[earlytype==1] + merr_log_n[earlytype==1]),
        	mag_sph[earlytype==1],0.5*(merr_mag_sph[earlytype==1] + perr_mag_sph[earlytype==1]),log_n[earlytype==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), mag_sph[earlytype==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), mag_sph[earlytype==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), mag_sph[earlytype==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), mag_sph[earlytype==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[earlytype==1])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['early'] == True:
		B_mfitexy = -24.74
		Berr_mfitexy = 0.14
		A_mfitexy = -8.99
		Aerr_mfitexy = 0.48
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[earlytype==1])),yy, color='red', ls='--', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[earlytype==1])),yy, color='red', ls='--', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[earlytype==1])), yy_lo, yy_up, alpha=0.2, facecolor='red')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[earlytype==1])), yy_lo, yy_up, alpha=0.2, facecolor='red')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_early.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
        	0.5*(perr_log_n[earlytype==1] + merr_log_n[earlytype==1]),
        	mag_sph[earlytype==1],0.5*(merr_mag_sph[earlytype==1] + perr_mag_sph[earlytype==1])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
        
	#######################################
	
        print 'bulges'
        print 'n', len(log_n[bulges==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[bulges==1]-np.average(log_n[bulges==1]),
        	0.5*(perr_log_n[bulges==1] + merr_log_n[bulges==1]),
        	mag_sph[bulges==1],0.5*(merr_mag_sph[bulges==1] + perr_mag_sph[bulges==1]),log_n[bulges==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[bulges==1]-np.average(log_n[bulges==1]), mag_sph[bulges==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[bulges==1]-np.average(log_n[bulges==1]), mag_sph[bulges==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[bulges==1]-np.average(log_n[bulges==1]), mag_sph[bulges==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[bulges==1]-np.average(log_n[bulges==1]), mag_sph[bulges==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[bulges==1])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['bulge'] == True:
		B_mfitexy = -22.18
		Berr_mfitexy = 0.20
		A_mfitexy = -4.34
		Aerr_mfitexy = 0.84
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[bulges==1])),yy, color='green', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[bulges==1])), yy_lo, yy_up, alpha=0.2, facecolor='green')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_bulge.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[bulges==1]-np.average(log_n[bulges==1]),
        	0.5*(perr_log_n[bulges==1] + merr_log_n[bulges==1]),
        	mag_sph[bulges==1],0.5*(merr_mag_sph[bulges==1] + perr_mag_sph[bulges==1])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
        
	#######################################
	
        print 'E'
        print 'n', len(log_n[simplemorphtype=='E'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']),
        	0.5*(perr_log_n[simplemorphtype=='E'] + merr_log_n[simplemorphtype=='E']),
        	mag_sph[simplemorphtype=='E'],0.5*(merr_mag_sph[simplemorphtype=='E'] + perr_mag_sph[simplemorphtype=='E']),log_n[simplemorphtype=='E']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), mag_sph[simplemorphtype=='E'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), mag_sph[simplemorphtype=='E'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), mag_sph[simplemorphtype=='E'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), mag_sph[simplemorphtype=='E'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='E'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['E'] == True:
		B_mfitexy = -25.74
		Berr_mfitexy = 0.19
		A_mfitexy = -10.07
		Aerr_mfitexy = 1.19
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='E'])),yy, color='red', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='E'])), yy_lo, yy_up, alpha=0.2, facecolor='red')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_ell.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']),
        	0.5*(perr_log_n[simplemorphtype=='E'] + merr_log_n[simplemorphtype=='E']),
        	mag_sph[simplemorphtype=='E'],0.5*(merr_mag_sph[simplemorphtype=='E'] + perr_mag_sph[simplemorphtype=='E'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
        
	#######################################
	
	print 'S0s'
        print 'n', len(log_n[simplemorphtype=='S0'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']),
        	0.5*(perr_log_n[simplemorphtype=='S0'] + merr_log_n[simplemorphtype=='S0']),
        	mag_sph[simplemorphtype=='S0'],0.5*(merr_mag_sph[simplemorphtype=='S0'] + perr_mag_sph[simplemorphtype=='S0']),log_n[simplemorphtype=='S0']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), mag_sph[simplemorphtype=='S0'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), mag_sph[simplemorphtype=='S0'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), mag_sph[simplemorphtype=='S0'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), mag_sph[simplemorphtype=='S0'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='S0'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['S0'] == True:
		B_mfitexy = -22.05
		Berr_mfitexy = 0.35
		A_mfitexy = -8.55
		Aerr_mfitexy = 2.79
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='S0'])),yy, color='darkorange', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='S0'])), yy_lo, yy_up, alpha=0.2, facecolor='darkorange')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_S0.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']),
        	0.5*(perr_log_n[simplemorphtype=='S0'] + merr_log_n[simplemorphtype=='S0']),
        	mag_sph[simplemorphtype=='S0'],0.5*(merr_mag_sph[simplemorphtype=='S0'] + perr_mag_sph[simplemorphtype=='S0'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
        
	#######################################
	
	print 'Sp'
        print 'n', len(log_n[simplemorphtype=='Sp'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']),
        	0.5*(perr_log_n[simplemorphtype=='Sp'] + merr_log_n[simplemorphtype=='Sp']),
        	mag_sph[simplemorphtype=='Sp'],0.5*(merr_mag_sph[simplemorphtype=='Sp'] + perr_mag_sph[simplemorphtype=='Sp']),log_n[simplemorphtype=='Sp']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), mag_sph[simplemorphtype=='Sp'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), mag_sph[simplemorphtype=='Sp'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), mag_sph[simplemorphtype=='Sp'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), mag_sph[simplemorphtype=='Sp'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='Sp'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['Sp'] == True:
		B_mfitexy = -22.23
		Berr_mfitexy = 0.33
		A_mfitexy = -3.60
		Aerr_mfitexy = 1.29
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])),yy, color='blue', ls='-', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])),yy, color='blue', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])), yy_lo, yy_up, alpha=0.2, facecolor='blue')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])), yy_lo, yy_up, alpha=0.2, facecolor='blue')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mag_sph_vs_n_maj_sph_spi.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn MAGsph err_MAGsph \n')
        for logNmaj, err_logNmaj, MAGsph, err_MAGsph in zip(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']),
        	0.5*(perr_log_n[simplemorphtype=='Sp'] + merr_log_n[simplemorphtype=='Sp']),
        	mag_sph[simplemorphtype=='Sp'],0.5*(merr_mag_sph[simplemorphtype=='Sp'] + perr_mag_sph[simplemorphtype=='Sp'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(MAGsph) + ' ' + str(err_MAGsph) + ' ' + '\n')
        datfile.close()
	
	##############################
        
	### akaike test
	X = log_n[all==1]-np.average(log_n[all==1])
	Y = mag_sph[all==1]
	err_Y = 0.5*(merr_mag_sph[all==1] + perr_mag_sph[all==1])
	m = -7.04
	q = -23.90
	k = 2
	print 'AICc(Y|X) single dataset:', akaike.get_AICc_singledataset(X,Y,err_Y,m,q,k)
	
	X1 = log_n[earlytype==1]-np.average(log_n[earlytype==1])
	Y1 = mag_sph[earlytype==1]
	err_Y1 = 0.5*(merr_mag_sph[earlytype==1] + perr_mag_sph[earlytype==1])
	m1 = -8.99
	q1 = -24.74
	X2 = log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp'])
	Y2 = mag_sph[simplemorphtype=='Sp']
	err_Y2 = 0.5*(merr_mag_sph[simplemorphtype=='Sp'] + perr_mag_sph[simplemorphtype=='Sp'])
	m2 = -3.60
	q2 = -22.23
	k = 4
	print 'AICc(Y|X) double dataset:', akaike.get_AICc_doubledataset(X1,Y1,err_Y1,m1,q1,X2,Y2,err_Y2,m2,q2,k)	
        
	##############################
        
	
	axislimits = [(0.3),(16),-19.01,-28.5]
	
	ax.set_xscale('log')
        ax.axis(axislimits)
        xticks = (np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
        
	ax_eb.set_xscale('log')
        ax_eb.axis(axislimits)
        ax_eb.set_xticks(xticks)
        ax_eb.set_xticklabels(xticks_labels)
        
	ax_eb.set_xlabel(r'$n_{\rm sph}^{\rm maj}$', labelpad=13)
        ax.set_ylabel(r'$MAG_{\rm sph}~[\rm mag]$', labelpad=13)
        ax_eb.set_ylabel(r'$MAG_{\rm sph}~[\rm mag]$', labelpad=13)
        plt.subplots_adjust(left=0.17,bottom=0.17,right=0.6,top=0.9)
        
	# legend
        markers_linlog.elliptical(ax, 'red', np.log10(0.4), (-27.9), 0.035, 0.12)
        ax.text(0.52, -27.6, 'E', fontsize=14)
        markers_linlog.lenticular(ax, 'red', np.log10(0.4), (-27.1), 0.03, 0.12)
        ax.text(0.52, -26.8, 'E/S0', fontsize=14)
        markers_linlog.lenticular(ax, 'darkorange', np.log10(0.4), (-26.3), 0.03, 0.12)
        ax.text(0.52, -26, 'S0', fontsize=14)
        markers_linlog.spiral(ax, 'darkorange', np.log10(1.2), (-27.9), 0.02, 0.07)
        ax.text(1.5, -27.6, 'S0/Sp', fontsize=14)
        markers_linlog.spiral(ax, 'blue', np.log10(1.2), (-27.1), 0.02, 0.07)
        ax.text(1.5, -26.8, 'Sp', fontsize=14)
        ax.scatter(1.2, -26.3, marker=r'$\star$', s=200, color='k', **scatter_kwargs) 
        ax.text(1.5, -26, 'merger', fontsize=14)

        # make inset
        #ins = plt.axes([.62, .18, .33, .27])
        #ins.axis([(0.3),(16),-19.01,-28.5])
	#ins.set_xscale('log')
        #ins.axes.get_xaxis().set_ticklabels([])
        #ins.axes.get_yaxis().set_ticklabels([])
        ax_eb.errorbar(n[simplemorphtype=='E'], mag_sph[simplemorphtype=='E'],
        	xerr=[merr_n[simplemorphtype=='E'],perr_n[simplemorphtype=='E']],
        	yerr=[merr_mag_sph[simplemorphtype=='E'],perr_mag_sph[simplemorphtype=='E']],
        	ecolor='red', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='E/S0'], mag_sph[simplemorphtype=='E/S0'],
        	xerr=[merr_n[simplemorphtype=='E/S0'],perr_n[simplemorphtype=='E/S0']],
        	yerr=[merr_mag_sph[simplemorphtype=='E/S0'],perr_mag_sph[simplemorphtype=='E/S0']],
        	ecolor='red', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='S0'], mag_sph[simplemorphtype=='S0'],
        	xerr=[merr_n[simplemorphtype=='S0'],perr_n[simplemorphtype=='S0']],
        	yerr=[merr_mag_sph[simplemorphtype=='S0'],perr_mag_sph[simplemorphtype=='S0']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='S0/Sp'], mag_sph[simplemorphtype=='S0/Sp'],
        	xerr=[merr_n[simplemorphtype=='S0/Sp'],perr_n[simplemorphtype=='S0/Sp']],
        	yerr=[merr_mag_sph[simplemorphtype=='S0/Sp'],perr_mag_sph[simplemorphtype=='S0/Sp']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='Sp'], mag_sph[simplemorphtype=='Sp'],
        	xerr=[merr_n[simplemorphtype=='Sp'],perr_n[simplemorphtype=='Sp']],
        	yerr=[merr_mag_sph[simplemorphtype=='Sp'],perr_mag_sph[simplemorphtype=='Sp']],
        	ecolor='blue', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='merger'], mag_sph[simplemorphtype=='merger'],
        	xerr=[merr_n[simplemorphtype=='merger'],perr_n[simplemorphtype=='merger']],
        	yerr=[merr_mag_sph[simplemorphtype=='merger'],perr_mag_sph[simplemorphtype=='merger']],
        	ecolor='k', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='out'], mag_sph[simplemorphtype=='out'],
        	xerr=[merr_n[simplemorphtype=='out'],perr_n[simplemorphtype=='out']],
        	yerr=[merr_mag_sph[simplemorphtype=='out'],perr_mag_sph[simplemorphtype=='out']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
         
	fig.subplots_adjust(hspace=0)
	
	#plt.show()
	#plt.savefig('test.pdf', format='pdf', dpi=100)
	plt.savefig(path_paper_figures + 'mag_vs_n_maj.pdf', format='pdf', dpi=1000)


def mbh_vs_n_sph_maj():
	
	outliers = [u'n0524', u'n3998']
	#outliers = []
	
	drawfit = {'all' : True, 'E' : False, 'S0' : False, 'Sp' : True, 'bulge' : False, 'early' : True}
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, \
		physres.log_n_maj_moffat_comb, \
		errV.perr_log_n, errV.merr_log_n, \
		anc.core \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	simplemorphtype = data[1]
	mbh = data[2].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[3].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[4].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	log_n = data[5].astype(np.float)
	n = 10**log_n
	perr_log_n = data[6].astype(np.float)
	merr_log_n = data[7].astype(np.float)
	perr_n = (10**perr_log_n - 1)*n
	merr_n = (10**(-perr_log_n) - 1)*(-n)
	core = data[8].astype(np.int)
		
        simplemorphtype[gal_id=='n0524'] = 'out'
	simplemorphtype[gal_id=='n3998'] = 'out'
	
	earlytype = mbh*[0]
        earlytype[simplemorphtype=='E'] = 1
        earlytype[simplemorphtype=='S0'] = 1
        earlytype[simplemorphtype=='E/S0'] = 1
	
	all = mbh*[0] + 1
	#all[simplemorphtype=='merger'] = 0
        all[gal_id=='n0524'] = 0
	all[gal_id=='n3998'] = 0
	
	# build figure
		       
	fig, (ax, ax_eb) = plt.subplots(2, sharex=True, sharey=True)

	scatter_kwargs = {"zorder":200}
	error_kwargs = {"lw":.5, "zorder":0}

	for x0,y0 in zip(n[simplemorphtype=='E'], mbh[simplemorphtype=='E']):
		markers_n.elliptical(ax, 'red', np.log10(x0), np.log10(y0), 0.032, 0.08)
	
 	for x0,y0 in zip(n[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
		markers_n.lenticular(ax, 'red', np.log10(x0), np.log10(y0), 0.03, 0.08)
       
	for x0,y0 in zip(n[simplemorphtype=='S0'], mbh[simplemorphtype=='S0']):
		markers_n.lenticular(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.03, 0.08)
	
        for x0,y0 in zip(n[simplemorphtype=='S0/Sp'], mbh[simplemorphtype=='S0/Sp']):
		markers_n.spiral(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.015, 0.04)
			
        for x0,y0 in zip(n[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
		markers_n.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.015, 0.04)
		
        #for x0,y0 in zip(n[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
	#	markers_n.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.03, 0.04)

	ax.scatter(n[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker=r'$\star$', s=200, color='k', **scatter_kwargs)	
			
        ax.text(n[gal_id=='n0524']-0.5, mbh[gal_id=='n0524']-1*10**8, 'N0524', color='darkorange', fontsize=12)
        ax.text(n[gal_id=='n3998']+0.15, mbh[gal_id=='n3998']-1*10**8, 'N3998', color='darkorange', fontsize=12)
	
	ax.scatter(n[gal_id=='n0524'], mbh[gal_id=='n0524'], marker='x', color='darkorange', s=100, linewidth=3)
	ax.scatter(n[gal_id=='n3998'], mbh[gal_id=='n3998'], marker='x', color='darkorange', s=100, linewidth=3)
	
	#expected_log_mbh = 8.18 + 3.39*(log_n - 0.51)
	#diff = log_mbh - expected_log_mbh 
	#twosout = gal_id[diff>2*0.58]
	#threesout = gal_id[diff>3*0.58]
	#print twosout, threesout

	#################################################
       
        print 'all'
        print 'n', len(log_n[all==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[all==1]-np.average(log_n[all==1]),
        	0.5*(perr_log_n[all==1] + merr_log_n[all==1]),
        	log_mbh[all==1],0.5*(merr_log_mbh[all==1] + perr_log_mbh[all==1]),log_n[all==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), log_mbh[all==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), log_mbh[all==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), log_mbh[all==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[all==1]-np.average(log_n[all==1]), log_mbh[all==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[all==1])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['all'] == True:
		B_mfitexy = 8.15 
		Berr_mfitexy = 0.06
		A_mfitexy = 3.37 
		Aerr_mfitexy = 0.15
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[all==1])),10**yy, color='k', ls='dotted', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[all==1])),10**yy, color='k', ls='dotted', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[all==1])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='k')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[all==1])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='k')
              
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mbh_vs_n_maj_sph_all.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn logMassBH err_logMassBH \n')
        for logNmaj, err_logNmaj, logMassBH, err_logMassBH in zip(log_n[all==1]-np.average(log_n[all==1]),
        	0.5*(perr_log_n[all==1] + merr_log_n[all==1]),
        	log_mbh[all==1],0.5*(merr_log_mbh[all==1] + perr_log_mbh[all==1])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
        datfile.close()
        
	#################################################
        
        print 'early'
        print 'n', len(log_n[earlytype==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
        	0.5*(perr_log_n[earlytype==1] + merr_log_n[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),log_n[earlytype==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), log_mbh[earlytype==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), log_mbh[earlytype==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), log_mbh[earlytype==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[earlytype==1]-np.average(log_n[earlytype==1]), log_mbh[earlytype==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[earlytype==1])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['early'] == True:
		B_mfitexy = 8.59 
		Berr_mfitexy = 0.07
		A_mfitexy = 3.58 
		Aerr_mfitexy = 0.27
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[earlytype==1])),10**yy, color='red', ls='--', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[earlytype==1])),10**yy, color='red', ls='--', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='red')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='red')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mbh_vs_n_maj_sph_early.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn logMassBH err_logMassBH \n')
        for logNmaj, err_logNmaj, logMassBH, err_logMassBH in zip(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
        	0.5*(perr_log_n[earlytype==1] + merr_log_n[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
        datfile.close()

        #################################################
        
        print 'Sp'
        print 'n', len(log_n[simplemorphtype=='Sp'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']),
        	0.5*(perr_log_n[simplemorphtype=='Sp'] + merr_log_n[simplemorphtype=='Sp']),
        	log_mbh[simplemorphtype=='Sp'],0.5*(merr_log_mbh[simplemorphtype=='Sp'] + perr_log_mbh[simplemorphtype=='Sp']),log_n[simplemorphtype=='Sp']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), log_mbh[simplemorphtype=='Sp'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), log_mbh[simplemorphtype=='Sp'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), log_mbh[simplemorphtype=='Sp'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']), log_mbh[simplemorphtype=='Sp'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='Sp'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['Sp'] == True:
		B_mfitexy = 7.24 
		Berr_mfitexy = 0.14
		A_mfitexy = 4.55 
		Aerr_mfitexy = 0.66
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])),10**yy, color='blue', ls='-', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])),10**yy, color='blue', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='blue')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='blue')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mbh_vs_n_maj_sph_spi.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn logMassBH err_logMassBH \n')
        for logNmaj, err_logNmaj, logMassBH, err_logMassBH in zip(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']),
        	0.5*(perr_log_n[simplemorphtype=='Sp'] + merr_log_n[simplemorphtype=='Sp']),
        	log_mbh[simplemorphtype=='Sp'],0.5*(merr_log_mbh[simplemorphtype=='Sp'] + perr_log_mbh[simplemorphtype=='Sp'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
        datfile.close()

        #################################################
        
        print 'S0'
        print 'n', len(log_n[simplemorphtype=='S0'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']),
        	0.5*(perr_log_n[simplemorphtype=='S0'] + merr_log_n[simplemorphtype=='S0']),
        	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']),log_n[simplemorphtype=='S0']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='S0'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['S0'] == True:
		B_mfitexy = 7.65  
		Berr_mfitexy = 0.12
		A_mfitexy = 3.78 
		Aerr_mfitexy = 0.85
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='S0'])),10**yy, color='darkorange', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='S0'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='darkorange')
              
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mbh_vs_n_maj_sph_S0.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn logMassBH err_logMassBH \n')
        for logNmaj, err_logNmaj, logMassBH, err_logMassBH in zip(log_n[simplemorphtype=='S0']-np.average(log_n[simplemorphtype=='S0']),
        	0.5*(perr_log_n[simplemorphtype=='S0'] + merr_log_n[simplemorphtype=='S0']),
        	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
        datfile.close()

        #################################################
        
        print 'E'
        print 'n', len(log_n[simplemorphtype=='E'])
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']),
        	0.5*(perr_log_n[simplemorphtype=='E'] + merr_log_n[simplemorphtype=='E']),
        	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E']),log_n[simplemorphtype=='E']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='E'])
        print
        print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0]), 'Delta =', "{0:.2f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1]), 'Delta =', "{0:.2f}".format(absscat_1)
        print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '    B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2]), 'Delta =', "{0:.2f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3]), 'Delta =', "{0:.2f}".format(absscat_3)
        print '---------------------------------'
       
        if drawfit['E'] == True:
		B_mfitexy = 8.91 
		Berr_mfitexy = 0.13
		A_mfitexy = 5.42 
		Aerr_mfitexy = 0.85
	
		logxx = np.arange(-10,20,0.1)
        	yy = (A_mfitexy*(logxx) + B_mfitexy)
        	ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='E'])),10**yy, color='red', ls='-', linewidth=2.)
        	ax_eb.plot(10**(logxx+np.average(log_n[simplemorphtype=='E'])),10**yy, color='red', ls='-', linewidth=2.)
        	#colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       
        	##### calculates 1sigma uncertainty band
        	yy_1 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_2 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy+Berr_mfitexy))
        	yy_3 = ((A_mfitexy+Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
        	yy_4 = ((A_mfitexy-Aerr_mfitexy)*(logxx) + (B_mfitexy-Berr_mfitexy))
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
        				
        	ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='E'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='red')
        	ax_eb.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='E'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='red')
       
        ### produce .dat file
        datfileName = '/Users/gsavorgnan/galaxy_vivisection/results/mbh_vs_n_maj_sph_ell.dat'
        datfile = open(datfileName, 'w')
        datfile.write('# lognmaj err_logn logMassBH err_logMassBH \n')
        for logNmaj, err_logNmaj, logMassBH, err_logMassBH in zip(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']),
        	0.5*(perr_log_n[simplemorphtype=='E'] + merr_log_n[simplemorphtype=='E']),
        	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E'])):
        	datfile.write(str(logNmaj) + ' ' + str(err_logNmaj) + ' ' + str(logMassBH) + ' ' + str(err_logMassBH) + ' ' + '\n')
        datfile.close()

        #################################################
        
	### akaike test
	X = log_n[all==1]-np.average(log_n[all==1])
	Y = log_mbh[all==1]
	err_Y = 0.5*(merr_log_mbh[all==1] + perr_log_mbh[all==1])
	m = 3.37 
	q = 8.15 
	k = 2
	print 'AICc(Y|X) single dataset:', akaike.get_AICc_singledataset(X,Y,err_Y,m,q,k)
	
	X1 = log_n[earlytype==1]-np.average(log_n[earlytype==1])
	Y1 = log_mbh[earlytype==1]
	err_Y1 = 0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1])
	m1 = 3.58 
	q1 = 8.59
	X2 = log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp'])
	Y2 = log_mbh[simplemorphtype=='Sp']
	err_Y2 = 0.5*(merr_log_mbh[simplemorphtype=='Sp'] + perr_log_mbh[simplemorphtype=='Sp'])
	m2 = 4.55
	q2 = 7.24 
	k = 4
	print 'AICc(Y|X) double dataset:', akaike.get_AICc_doubledataset(X1,Y1,err_Y1,m1,q1,X2,Y2,err_Y2,m2,q2,k)	
        
	##############################

	# legend
        markers_n.elliptical(ax, 'red', np.log10(0.4), 11.00, 0.032, 0.08)
        ax.text(0.52, 10**10.85, 'E', fontsize=14)
        markers_n.lenticular(ax, 'red', np.log10(0.4), 10.50, 0.03, 0.08)
        ax.text(0.52, 10**10.35, 'E/S0', fontsize=14)
        markers_n.lenticular(ax, 'darkorange', np.log10(0.4), 10.00, 0.03, 0.08)
        ax.text(0.52, 10**9.85, 'S0', fontsize=14)
        markers_n.spiral(ax, 'darkorange', np.log10(1.2), 11.00, 0.015, 0.04)
        ax.text(1.45, 10**10.85, 'S0/Sp', fontsize=14)
        markers_n.spiral(ax, 'blue', np.log10(1.2), 10.50, 0.015, 0.04)
        ax.text(1.45, 10**10.35, 'Sp', fontsize=14)
        ax.scatter([1.2], [10**10.00], marker=r'$\star$', s=200, color='k', **scatter_kwargs) 
        ax.text(1.45, 10**9.85, 'merger', fontsize=14)

	axislimits = [(0.3),(16),10**5.05,10**11.5]
	
	ax.set_xscale('log')
	ax.set_yscale('log')
        ax.axis(axislimits)
        xticks = (np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
        
	ax_eb.set_xscale('log')
	ax_eb.set_yscale('log')
        ax_eb.axis(axislimits)
        ax_eb.set_xticks(xticks)
        ax_eb.set_xticklabels(xticks_labels)
        
	
	ax_eb.set_xlabel(r'$n_{\rm sph}^{\rm maj}$', labelpad=13)
        ax.set_ylabel(r'$M_{\rm BH}~[\rm M_\odot]$', labelpad=13)
        ax_eb.set_ylabel(r'$M_{\rm BH}~[\rm M_\odot]$', labelpad=13)
	plt.subplots_adjust(left=0.15,bottom=0.17,right=0.6,top=0.9)
        
        # make inset
       #ins = plt.axes([.61, .18, .34, .3])
       #ins.axis([(0.3),(16),10**5.5,10**10.8])
       #ins.set_xscale('log')
       #ins.set_yscale('log')
       #ins.axes.get_xaxis().set_ticklabels([])
       #ins.axes.get_yaxis().set_ticklabels([])
        ax_eb.errorbar(n[simplemorphtype=='E'], mbh[simplemorphtype=='E'],
        	xerr=[merr_n[simplemorphtype=='E'],perr_n[simplemorphtype=='E']],
        	yerr=[merr_mbh[simplemorphtype=='E'],perr_mbh[simplemorphtype=='E']],
        	ecolor='red', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'],
        	xerr=[merr_n[simplemorphtype=='E/S0'],perr_n[simplemorphtype=='E/S0']],
        	yerr=[merr_mbh[simplemorphtype=='E/S0'],perr_mbh[simplemorphtype=='E/S0']],
        	ecolor='red', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'],
        	xerr=[merr_n[simplemorphtype=='S0'],perr_n[simplemorphtype=='S0']],
        	yerr=[merr_mbh[simplemorphtype=='S0'],perr_mbh[simplemorphtype=='S0']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='S0/Sp'], mbh[simplemorphtype=='S0/Sp'],
        	xerr=[merr_n[simplemorphtype=='S0/Sp'],perr_n[simplemorphtype=='S0/Sp']],
        	yerr=[merr_mbh[simplemorphtype=='S0/Sp'],perr_mbh[simplemorphtype=='S0/Sp']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp'],
        	xerr=[merr_n[simplemorphtype=='Sp'],perr_n[simplemorphtype=='Sp']],
        	yerr=[merr_mbh[simplemorphtype=='Sp'],perr_mbh[simplemorphtype=='Sp']],
        	ecolor='blue', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'],
        	xerr=[merr_n[simplemorphtype=='merger'],perr_n[simplemorphtype=='merger']],
        	yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']],
        	ecolor='k', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
        ax_eb.errorbar(n[simplemorphtype=='out'], mbh[simplemorphtype=='out'],
        	xerr=[merr_n[simplemorphtype=='out'],perr_n[simplemorphtype=='out']],
        	yerr=[merr_mbh[simplemorphtype=='out'],perr_mbh[simplemorphtype=='out']],
        	ecolor='darkorange', marker='', ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

	fig.subplots_adjust(hspace=0)
       	
        #plt.show()
	#plt.savefig('test.pdf', format='pdf', dpi=100)
	plt.savefig(path_paper_figures + 'mbh_vs_n_maj.pdf', format='pdf', dpi=1000)

	
		
def main():
	mag_sph_vs_n_maj()
	#mag_sph_vs_n_eq()
	#mbh_vs_n_sph_maj()
	#mbh_vs_n_sph_eq()
	
main()		
		
	
