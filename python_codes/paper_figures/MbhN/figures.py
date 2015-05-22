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

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 26})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/MbhN/images/'



def mbh_vs_n_sph(axis):
	
	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, \
		physres.log_n_' + axis + '_moffat_comb, \
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
		
       #ELL_core = ELL*[0]
       #ELL_sersic = ELL*[0]
       #BUL_core = ELL*[0]
       #BUL_sersic = ELL*[0]
       #
        earlytype = mbh*[0]
        earlytype[simplemorphtype=='E'] = 1
        earlytype[simplemorphtype=='S0'] = 1
        earlytype[simplemorphtype=='E/S0'] = 1
       
       #bulge = core*[0]
       #bulge[simplemorphtype=='S0'] = 1
       #bulge[simplemorphtype=='S'] = 1
       #bulge[simplemorphtype=='S0/S'] = 1
       #
       #for i in range(len(ELL)):
       #	if simplemorphtype[i]=='E' and core[i]==1 and gal_id[i] not in outliers:
       #		ELL_core[i] = 1
       #	elif simplemorphtype[i]=='E' and core[i]==0 and gal_id[i] not in outliers:
       #		ELL_sersic[i] = 1
       #	elif bulge[i]==1 and core[i]==1 and gal_id[i] not in outliers:
       #		BUL_core[i] = 1
       #	elif bulge[i]==1 and core[i]==0 and gal_id[i] not in outliers:
       #		BUL_sersic[i] = 1
       #
       #morph_coreList = []
       #for i in range(len(simplemorphtype)):
       #	if gal_id[i] not in outliers:
       #		morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
       #	else:
       #		morph_coreList.append('outlier')	
       #morph_core = np.asarray(morph_coreList)
	       
        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":200}
	error_kwargs = {"lw":.5, "zorder":0}

	for x0,y0 in zip(n[simplemorphtype=='E'], mbh[simplemorphtype=='E']):
		markers_n.elliptical(ax, 'red', np.log10(x0), np.log10(y0), 0.032, 0.06)
	
 	for x0,y0 in zip(n[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
		markers_n.lenticular(ax, 'red', np.log10(x0), np.log10(y0), 0.03, 0.06)
       
	for x0,y0 in zip(n[simplemorphtype=='S0'], mbh[simplemorphtype=='S0']):
		markers_n.lenticular(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.03, 0.06)
	
        for x0,y0 in zip(n[simplemorphtype=='S0/Sp'], mbh[simplemorphtype=='S0/Sp']):
		markers_n.spiral(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.015, 0.03)
			
        for x0,y0 in zip(n[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
		markers_n.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.015, 0.03)
		
        #for x0,y0 in zip(n[simplemorphtype=='Sp'], mbh[simplemorphtype=='Sp']):
	#	markers_n.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.03, 0.04)

	ax.scatter(n[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
			
       ##################################################
       #
       #print 'all'
       #print 'n', len(log_n)
       #A,B,Aerr,Berr,covAB=bces.bces(log_n-np.average(log_n),
       #	0.5*(perr_log_n + merr_log_n),
       #	log_mbh,0.5*(merr_log_mbh + perr_log_mbh),log_n*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n)
       #print
       ##print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1])
       #print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_n)),10**yy, color='k', ls='-', linewidth=2.)
       ##colorline.colorline(10**(logxx+np.average(log_n)), 10**yy, cmap=green_red)
       #
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
       #ax.fill_between(10**(logxx+np.average(log_n)), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='k')
       #
       ##################################################
       #
       #print 'early'
       #print 'n', len(log_n[earlytype==1])
       #A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
       #	0.5*(perr_log_n[earlytype==1] + merr_log_n[earlytype==1]),
       #	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),log_n[earlytype==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n[earlytype==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1])
       #print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_n[earlytype==1])),10**yy, color='r', ls='--', linewidth=2.)
       ##colorline.colorline(10**(logxx+np.average(log_n[earlytype==1])), 10**yy, cmap=green_red)
       #
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
       #ax.fill_between(10**(logxx+np.average(log_n[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.3, facecolor='r')
       #
       ##################################################
       #
       #print 'late'
       #print 'n', len(log_n[simplemorphtype=='Sp'])
       #A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='Sp']-np.average(log_n[simplemorphtype=='Sp']),
       #	0.5*(perr_log_n[simplemorphtype=='Sp'] + merr_log_n[simplemorphtype=='Sp']),
       #	log_mbh[simplemorphtype=='Sp'],0.5*(merr_log_mbh[simplemorphtype=='Sp'] + perr_log_mbh[simplemorphtype=='Sp']),log_n[simplemorphtype=='Sp']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n[simplemorphtype=='Sp'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.2f}".format(A[0]), '\pm', "{0:.2f}".format(Aerr[0]), '	B = ', "{0:.2f}".format(B[0]), '\pm', "{0:.2f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.2f}".format(A[1]), '\pm', "{0:.2f}".format(Aerr[1]), '	B = ', "{0:.2f}".format(B[1]), '\pm', "{0:.2f}".format(Berr[1])
       #print 'bisector    A =', "{0:.2f}".format(A[2]), '\pm', "{0:.2f}".format(Aerr[2]), '   B = ', "{0:.2f}".format(B[2]), '\pm', "{0:.2f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.2f}".format(A[3]), '\pm', "{0:.2f}".format(Aerr[3]), '	B = ', "{0:.2f}".format(B[3]), '\pm', "{0:.2f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])),10**yy, color='b', ls='-', linewidth=2.)
       #
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
       #ax.fill_between(10**(logxx+np.average(log_n[simplemorphtype=='Sp'])), 10**yy_lo, 10**yy_up, alpha=0.2, facecolor='b')
       #
       ## legend
       #markers.elliptical(ax, 'red', 8.85, 10.9, 0.08)
       #ax.text(10**9.1, 10**10.75, 'E')
       #markers.lenticular(ax, 'red', 8.85, 10.45, 0.08)
       #ax.text(10**9.1, 10**10.3, 'E/S0')
       #markers.lenticular(ax, 'darkorange', 8.85, 10., 0.08)
       #ax.text(10**9.1, 10**9.85, 'S0')
       #markers.spiral(ax, 'darkorange', 9.8, 10.9, 0.04)
       #ax.text(10**10.05, 10**10.75, 'S0/Sp')
       #markers.spiral(ax, 'blue', 9.8, 10.45, 0.04)
       #ax.text(10**10.05, 10**10.3, 'Sp')
       #ax.scatter([10**9.8], [10**10.], marker=r'$\star$', s=500, color='k', **scatter_kwargs) 
       #ax.text(10**10.05, 10**9.85, 'merger')
	
	ax.set_xscale('log')
	ax.set_yscale('log')
        ax.axis([(0.3),(18),10**5.3,10**11.2])
        xticks = (np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
        plt.xlabel(r'spheroid central radial concentration', labelpad=13)
        plt.ylabel(r'black hole mass $[\rm M_\odot]$', labelpad=13)
	plt.subplots_adjust(left=0.15,bottom=0.17,right=0.97,top=0.9)
        #plt.show()
	plt.savefig(path_paper_figures + 'mbh_vs_n_' + axis + '.pdf', format='pdf', dpi=1000)

		
		
def main():
	mbh_vs_n_sph('maj')

main()		
		
	
