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
from instruments.linear_regression import akaike

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 22})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_auxiliary_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/auxiliary/'
path_scalrel_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/scaling_relations/'



def make_colormap(colors):
    """
    Define a new color map based on values specified in the dictionary
    colors, where colors[z] is the color that value z should be mapped to,
    with linear interpolation between the given values of z.

    The z values (dictionary keys) are real numbers and the values
    colors[z] can be either an RGB list, e.g. [1,0,0] for red, or an
    html hex string, e.g. "#ff0000" for red.
    """

    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort
    
    z = sort(colors.keys())
    n = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)
    
    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        #i'th color at level z[i]:
        Ci = colors[z[i]]      
        if type(Ci) == str:
            # a hex string of form '#ff0000' for example (for red)
            RGB = CC.to_rgb(Ci)
        else:
            # assume it's an RGB triple already:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])

    cmap_dict = {}
    cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap


red_green = make_colormap({0.:'r', 1.:'g'})

green_red = make_colormap({0.:'g', 1.:'r'})


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


def mag_sph_vs_n(axis):
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, \
		pysres.mag_sph_eq_moffat_comb, pysres.log_n_' + axis + '_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, errV.perr_log_n, errV.merr_log_n, \
		pysres.mag_tot_eq_moffat_comb \
        	FROM Ancillary AS anc \
        	JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
        	WHERE anc.fit1D_done = 1;'
       
        cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
        #print data
        simplemorphtype = data[1]
        mag_sph = data[2].astype(np.float)
        log_n = data[3].astype(np.float)
	perr_mag_sph = data[4].astype(np.float)
	merr_mag_sph = data[5].astype(np.float)
	perr_log_n = data[6].astype(np.float)
	merr_log_n = data[7].astype(np.float)
	mag_tot = data[8].astype(np.float)
	
	BT = 10**(-0.4*(mag_sph))/10**(-0.4*(mag_tot))
	
	earlytype = mag_sph*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = mag_sph*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
       
        fig, ax = plt.subplots()
       #ax.errorbar(log_n[simplemorphtype=='E'], mag_sph[simplemorphtype=='E'], 
       #	xerr=[merr_log_n[simplemorphtype=='E'],perr_log_n[simplemorphtype=='E']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='E'],perr_mag_sph[simplemorphtype=='E']], 
       #	ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
       #ax.errorbar(log_n[simplemorphtype=='S0'], mag_sph[simplemorphtype=='S0'], 
       #	xerr=[merr_log_n[simplemorphtype=='S0'],perr_log_n[simplemorphtype=='S0']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='S0'],perr_mag_sph[simplemorphtype=='S0']], 
       #	ecolor='green', fmt='go', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
       #ax.errorbar(log_n[simplemorphtype=='E/S0'], mag_sph[simplemorphtype=='E/S0'], 
       #	xerr=[merr_log_n[simplemorphtype=='E/S0'],perr_log_n[simplemorphtype=='E/S0']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='E/S0'],perr_mag_sph[simplemorphtype=='E/S0']], 
       #	ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
       #ax.errorbar(log_n[simplemorphtype=='S'], mag_sph[simplemorphtype=='S'], 
       #	xerr=[merr_log_n[simplemorphtype=='S'],perr_log_n[simplemorphtype=='S']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='S'],perr_mag_sph[simplemorphtype=='S']], 
       #	ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
       #ax.errorbar(log_n[simplemorphtype=='S0/S'], mag_sph[simplemorphtype=='S0/S'], 
       #	xerr=[merr_log_n[simplemorphtype=='S0/S'],perr_log_n[simplemorphtype=='S0/S']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='S0/S'],perr_mag_sph[simplemorphtype=='S0/S']], 
       #	ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
       #ax.errorbar(log_n[simplemorphtype=='merger'], mag_sph[simplemorphtype=='merger'], 
       #	xerr=[merr_log_n[simplemorphtype=='merger'],perr_log_n[simplemorphtype=='merger']], 
       #	yerr=[merr_mag_sph[simplemorphtype=='merger'],perr_mag_sph[simplemorphtype=='merger']], 
       #	ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
       
        ax.scatter(log_n, mag_sph, c=BT, s=80)

        xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
	
	### fit using BCES ####
        print 'early'
        #print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
        A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),0.5*(perr_log_n[earlytype==1]+merr_log_n[earlytype==1]),mag_sph[earlytype==1],
        	0.5*(merr_mag_sph[earlytype==1] + perr_mag_sph[earlytype==1]),log_n[earlytype==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[earlytype==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
        
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        #print A[2], B[2]
        # calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_n[ELL==1]-np.average(log_n[ELL==1]),mag_sph[ELL==1],A[2],B[2],conf=0.68,x=logxx)
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
       
        print 'ell'
	print log_n[simplemorphtype=='E'], mag_sph[simplemorphtype=='E']
        #print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']),0.5*(perr_log_n[simplemorphtype=='E']+merr_log_n[simplemorphtype=='E']),mag_sph[simplemorphtype=='E'],
        	0.5*(merr_mag_sph[simplemorphtype=='E'] + perr_mag_sph[simplemorphtype=='E']),log_n[simplemorphtype=='E']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='E'])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
        
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        #print A[2], B[2]
        # calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_n[ELL==1]-np.average(log_n[ELL==1]),mag_sph[ELL==1],A[2],B[2],conf=0.68,x=logxx)
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
       
        #ax.plot(logxx+np.average(log_n[simplemorphtype=='E']),yy, color='r', ls='-', linewidth=2.)
        # plots a shaded area containing the prediction band  
        #plt.fill_between(logxx+np.average(log_n[ELL==1]), lpb68, upb68, alpha=0.15, facecolor='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),lpb68, color='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),upb68, color='red')
       
	### fit using BCES ####
        print 'bulge'
        #print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
        A,B,Aerr,Berr,covAB=bces.bces(log_n[bulge==1]-np.average(log_n[bulge==1]),0.5*(perr_log_n[bulge==1]+merr_log_n[bulge==1]),mag_sph[bulge==1],
        	0.5*(merr_mag_sph[bulge==1] + perr_mag_sph[bulge==1]),log_n[bulge==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[bulge==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
        
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        #print A[2], B[2]
        # calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_n[ELL==1]-np.average(log_n[ELL==1]),mag_sph[ELL==1],A[2],B[2],conf=0.68,x=logxx)
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
       
        #ax.plot(logxx+np.average(log_n[bulge==1]),yy, color='y', ls='-', linewidth=2.)
        # plots a shaded area containing the prediction band  
        #plt.fill_between(logxx+np.average(log_n[ELL==1]), lpb68, upb68, alpha=0.15, facecolor='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),lpb68, color='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),upb68, color='red')
       
        print 'spi'
        #print log_n-np.average(log_n),0.5*(perr_log_n+merr_log_n),mag_sph, 0.5*(merr_mag_sph + perr_mag_sph),log_n*[0.0]
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='S']-np.average(log_n[simplemorphtype=='S']),0.5*(perr_log_n[simplemorphtype=='S']+merr_log_n[simplemorphtype=='S']),mag_sph[simplemorphtype=='S'],
        	0.5*(merr_mag_sph[simplemorphtype=='S'] + perr_mag_sph[simplemorphtype=='S']),log_n[simplemorphtype=='S']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='S'])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
        
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        #print A[2], B[2]
        # calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_n[ELL==0]-np.average(log_n[ELL==0]),mag_sph[ELL==0],A[2],B[2],conf=0.68,x=logxx)
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
       
        #ax.plot(logxx+np.average(log_n[simplemorphtype=='S']),yy, color='b', ls='-', linewidth=2.)
        # plots a shaded area containing the prediction band  
        #plt.fill_between(logxx+np.average(log_n[ELL==1]), lpb68, upb68, alpha=0.15, facecolor='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),lpb68, color='red')
        #ax.plot(logxx+np.average(log_n[ELL==1]),upb68, color='red')
       
        ### fit using FITEXY ###
        print 'Ellipticals'
        #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_n[simplemorphtype=='E']-np.average(log_n[simplemorphtype=='E']), 0.5*(perr_log_n[simplemorphtype=='E']+merr_log_n[simplemorphtype=='E']),
        #	mag_sph[simplemorphtype=='E'], 0.5*(merr_mag_sph[simplemorphtype=='E'] + perr_mag_sph[simplemorphtype=='E']))
        xx = np.arange(0.001,30,0.1)
        logxx = np.log10(xx)
        # plot bisector relation
        #y_bisec = A + B*logxx
        y_bisec = -25.71 -7.39*logxx
        #ax.plot(logxx+np.average(log_n[simplemorphtype=='E']), y_bisec, ls='--', color='red', linewidth=2.)
       
       #### fit using FITEXY ###
       #print 'Bulges'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_n[ELL==0]-np.average(log_n[ELL==0]), 0.5*(perr_log_n[ELL==0]+merr_log_n[ELL==0]),
       #	mag_sph[ELL==0], 0.5*(merr_mag_sph[ELL==0] + perr_mag_sph[ELL==0]))
       #xx = np.arange(0.001,30,0.1)
       #logxx = np.log10(xx)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(log_n[ELL==0]), y_bisec, ls='-', color='blue', linewidth=2.)
	
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
        plt.xlabel(r'$n^{\rm ' + axis + '}$', labelpad=15)
        plt.ylabel(r'$MAG_{\rm sph} \rm ~[mag]$', labelpad=15)
        plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
        #plt.savefig(path_scalrel_plots + 'mag_sph_vs_n' + axis + '.pdf', format='pdf', dpi=1000)
		
	
	

def mbh_vs_n(axis):
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.simplemorphtype, anc.core, \
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
        simplemorphtype = data[4]
	core = data[5].astype(np.int)
	log_n = data[6].astype(np.float)
	n = 10**log_n
	perr_log_n = data[7].astype(np.float)
	merr_log_n = data[8].astype(np.float)
	
	earlytype = core*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = core*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
		
        morph_coreList = []
        for i in range(len(simplemorphtype)):
        	morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
        morph_core = np.asarray(morph_coreList)

        fig, ax = plt.subplots()
        plt.axis([np.log10(0.25),np.log10(18),10**5.5,10**10.8])
	
	print 'early'
        A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),0.5*(perr_log_n[earlytype==1]+merr_log_n[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),log_n[earlytype==1]*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[earlytype==1])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
	logxx = np.arange(np.log10(0.1),np.log10(5.2),0.01)
        yy = (A[2]*(logxx) + B[2])
        colorline.colorline(logxx+np.average(log_n[earlytype==1]), 10**yy, cmap=green_red)
	logxx = np.arange(np.log10(0.01),np.log10(5.2),0.01)
        ##### calculates the prediction bands for the given input arrays
        lpb68_g,upb68_g,logxx_g = predband.predband(log_n[earlytype==1]-np.average(log_n[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.68,x=logxx[logxx+np.average(log_n[earlytype==1])<0.5])
        lpb68_r,upb68_r,logxx_r = predband.predband(log_n[earlytype==1]-np.average(log_n[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.68,x=logxx[logxx+np.average(log_n[earlytype==1])>=0.49])
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        #### plots a shaded area containing the prediction band  
        ax.fill_between(logxx_g+np.average(log_n[earlytype==1]), 10**lpb68_g, 10**upb68_g, alpha=0.15, facecolor='green')
        ax.fill_between(logxx_r+np.average(log_n[earlytype==1]), 10**lpb68_r, 10**upb68_r, alpha=0.15, facecolor='red')
	
	print 'spi'
        A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='S']-np.average(log_n[simplemorphtype=='S']),0.5*(perr_log_n[simplemorphtype=='S']+merr_log_n[simplemorphtype=='S']),
        	log_mbh[simplemorphtype=='S'],0.5*(merr_log_mbh[simplemorphtype=='S'] + perr_log_mbh[simplemorphtype=='S']),log_n[simplemorphtype=='S']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[simplemorphtype=='S'])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
	logxx = np.arange(np.log10(0.1),np.log10(20.2),0.01)
        yy = (A[2]*(logxx) + B[2])
        
	ax.plot(logxx+np.average(log_n[simplemorphtype=='S']),10**yy, color='blue', ls='--', linewidth=4.)
	##### calculates the prediction bands for the given input arrays
        lpb68,upb68,logxx = predband.predband(log_n[simplemorphtype=='S']-np.average(log_n[simplemorphtype=='S']),log_mbh[simplemorphtype=='S'],A[2],B[2],conf=0.68,x=logxx)
        ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
        #### plots a shaded area containing the prediction band  
        ax.fill_between(logxx+np.average(log_n[simplemorphtype=='S']), 10**lpb68, 10**upb68, alpha=0.15, facecolor='blue')
	
			
        ax.errorbar(log_n[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_log_n[simplemorphtype=='merger'],perr_log_n[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)

	ax.errorbar(log_n[morph_core=='E_1'], mbh[morph_core=='E_1'], 
		xerr=[merr_log_n[morph_core=='E_1'],perr_log_n[morph_core=='E_1']], 
		yerr=[merr_mbh[morph_core=='E_1'],perr_mbh[morph_core=='E_1']], 
		ecolor='red', marker='o', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(log_n[morph_core==u'E_0'], mbh[morph_core==u'E_0'], 
		xerr=[merr_log_n[morph_core==u'E_0'],perr_log_n[morph_core==u'E_0']], 
		yerr=[merr_mbh[morph_core==u'E_0'],perr_mbh[morph_core==u'E_0']], 
		ecolor='red', marker='o', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(log_n[morph_core=='S0_1'], mbh[morph_core=='S0_1'], 
		xerr=[merr_log_n[morph_core=='S0_1'],perr_log_n[morph_core=='S0_1']], 
		yerr=[merr_mbh[morph_core=='S0_1'],perr_mbh[morph_core=='S0_1']], 
		ecolor='green', marker='o', mfc='white', mec='green', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(log_n[morph_core=='S0_0'], mbh[morph_core=='S0_0'], 
		xerr=[merr_log_n[morph_core=='S0_0'],perr_log_n[morph_core=='S0_0']], 
		yerr=[merr_mbh[morph_core=='S0_0'],perr_mbh[morph_core=='S0_0']], 
		ecolor='green', marker='o', mfc='green', mec='green', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(log_n[morph_core=='E/S0_1'], mbh[morph_core=='E/S0_1'], 
		xerr=[merr_log_n[morph_core=='E/S0_1'],perr_log_n[morph_core=='E/S0_1']], 
		yerr=[merr_mbh[morph_core=='E/S0_1'],perr_mbh[morph_core=='E/S0_1']], 
		ecolor='red', marker='*', mfc='white', mec='red', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(log_n[morph_core=='E/S0_0'], mbh[morph_core=='E/S0_0'], 
		xerr=[merr_log_n[morph_core=='E/S0_0'],perr_log_n[morph_core=='E/S0_0']], 
		yerr=[merr_mbh[morph_core=='E/S0_0'],perr_mbh[morph_core=='E/S0_0']], 
		ecolor='red', marker='*', mfc='red', mec='red', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(log_n[morph_core=='S_1'], mbh[morph_core=='S_1'], 
        	xerr=[merr_log_n[morph_core=='S_1'],perr_log_n[morph_core=='S_1']], 
        	yerr=[merr_mbh[morph_core=='S_1'],perr_mbh[morph_core=='S_1']], 
        	ecolor='blue', marker='o', mfc='white', mec='blue', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(log_n[morph_core=='S_0'], mbh[morph_core=='S_0'], 
        	xerr=[merr_log_n[morph_core=='S_0'],perr_log_n[morph_core=='S_0']], 
        	yerr=[merr_mbh[morph_core=='S_0'],perr_mbh[morph_core=='S_0']], 
        	ecolor='blue', marker='o', mfc='blue', mec='blue', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(log_n[morph_core=='S0/S_1'], mbh[morph_core=='S0/S_1'], 
		xerr=[merr_log_n[morph_core=='S0/S_1'],perr_log_n[morph_core=='S0/S_1']], 
		yerr=[merr_mbh[morph_core=='S0/S_1'],perr_mbh[morph_core=='S0/S_1']], 
		ecolor='blue', marker='*', mfc='white', mec='blue', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(log_n[morph_core=='S0/S_0'], mbh[morph_core=='S0/S_0'], 
		xerr=[merr_log_n[morph_core=='S0/S_0'],perr_log_n[morph_core=='S0/S_0']], 
		yerr=[merr_mbh[morph_core=='S0/S_0'],perr_mbh[morph_core=='S0/S_0']], 
		ecolor='blue', marker='*', mfc='blue', mec='blue', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)


        ax.set_yscale('log')
        xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
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
        plt.xlabel(r'$n^{\rm ' + axis + '}$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_n_' + axis + '.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mag_sph():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
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
	
        earlytype = ELL*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = core*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
	
        ELL_core = ELL*[0]
        ELL_sersic = ELL*[0]
        BUL_core = ELL*[0]
        BUL_sersic = ELL*[0]
       
        for i in range(len(ELL)):
        	if simplemorphtype[i]=='E' and core[i]==1:
        		ELL_core[i] = 1
        	elif simplemorphtype[i]=='E' and core[i]==0:
        		ELL_sersic[i] = 1
        	elif bulge[i]==1 and core[i]==1:
        		BUL_core[i] = 1
        	elif bulge[i]==1 and core[i]==0:
        		BUL_sersic[i] = 1
        
	X = mag_sph[bulge==1]
	Y = log_mbh[bulge==1]
	C = core[bulge==1]
	print len(X[C==0])
	print stats.spearmanr(X[C==0], Y[C==0])
	
        morph_coreList = []
        for i in range(len(simplemorphtype)):
        	morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
        morph_core = np.asarray(morph_coreList)
	
        fig, ax = plt.subplots()
        plt.axis([-19.01,-27.99,10**5.5,10**10.8])

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
       #print 'ser BUL'
       #print '<mag_sph>', np.average(mag_sph[BUL_sersic==1])
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[BUL_sersic==1]-np.average(mag_sph[BUL_sersic==1]), 0.5*(perr_mag_sph[BUL_sersic==1]+merr_mag_sph[BUL_sersic==1]),
       #	log_mbh[BUL_sersic==1], 0.5*(merr_log_mbh[BUL_sersic==1] + perr_log_mbh[BUL_sersic==1]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[BUL_sersic==1]), 10**y_bisec, ls='-', color='pink', linewidth=2.)
       ##label_bul = r'$B_{\rm bul} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-24, 10**6, label_bul, fontsize=20)
       
       #### fit using FITEXY ###
       #print 'spi'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[simplemorphtype=='S']-np.average(mag_sph[simplemorphtype=='S']), 0.5*(perr_mag_sph[simplemorphtype=='S']+merr_mag_sph[simplemorphtype=='S']),
       #	log_mbh[simplemorphtype=='S'], 0.5*(merr_log_mbh[simplemorphtype=='S'] + perr_log_mbh[simplemorphtype=='S']))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_sph[simplemorphtype=='S']), 10**y_bisec, ls='-', color='blue', linewidth=2.)
       ##label_bul = r'$B_{\rm bul} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-24, 10**6, label_bul, fontsize=20)
       
       #ax.errorbar(mag_sph[ELL_core==1], mbh[ELL_core==1], xerr=[merr_mag_sph[ELL_core==1],perr_mag_sph[ELL_core==1]], yerr=[merr_mbh[ELL_core==1],perr_mbh[ELL_core==1]], ecolor='red', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[ELL_sersic==1], mbh[ELL_sersic==1], xerr=[merr_mag_sph[ELL_sersic==1],perr_mag_sph[ELL_sersic==1]], yerr=[merr_mbh[ELL_sersic==1],perr_mbh[ELL_sersic==1]], ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[BUL_core==1], mbh[BUL_core==1], xerr=[merr_mag_sph[BUL_core==1],perr_mag_sph[BUL_core==1]], yerr=[merr_mbh[BUL_core==1],perr_mbh[BUL_core==1]], ecolor='blue', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
       #ax.errorbar(mag_sph[BUL_sersic==1], mbh[BUL_sersic==1], xerr=[merr_mag_sph[BUL_sersic==1],perr_mag_sph[BUL_sersic==1]], yerr=[merr_mbh[BUL_sersic==1],perr_mbh[BUL_sersic==1]], ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
      # print 'core'
      # A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==1]-np.average(mag_sph[core==1]),0.5*(perr_mag_sph[core==1]+merr_mag_sph[core==1]),
      # 	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),mag_sph[core==1]*[0.0])
      # print '---------------------------------'
      # print 'y = A*(x-<x>) + B '
      # print '<x> =', np.average(mag_sph[core==1])
      # print
      # #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      # #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      # print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      # #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      # print '---------------------------------'
      #
      # logxx = np.arange(-40,40,0.01)
      # yy = (A[2]*(logxx) + B[2])
      # ax.plot(logxx+np.average(mag_sph[core==1]),10**yy, color='k', ls='--', linewidth=2.)
      # m1 = A[2]
      # q1 = B[2]
      #
      # print 'Sersic'
      # A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==0]-np.average(mag_sph[core==0]),0.5*(perr_mag_sph[core==0]+merr_mag_sph[core==0]),
      # 	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),mag_sph[core==0]*[0.0])
      # print '---------------------------------'
      # print 'y = A*(x-<x>) + B '
      # print '<x> =', np.average(mag_sph[core==0])
      # print
      # #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
      # #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
      # print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
      # #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
      # print '---------------------------------'
      #
      # logxx = np.arange(-40,40,0.01)
      # yy = (A[2]*(logxx) + B[2])
      # ax.plot(logxx+np.average(mag_sph[core==0]),10**yy, color='k', ls='-', linewidth=2.)
      # m2 = A[2]
      # q2 = B[2]
      #
      # AICc_core = akaike.get_AICc_doubledataset(mag_sph[core==1]-np.average(mag_sph[core==1]),log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),m1,q1, 
      # 	mag_sph[core==0]-np.average(mag_sph[core==0]),log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),m2,q2,4)
       
        print 'all'
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph-np.average(mag_sph),0.5*(perr_mag_sph+merr_mag_sph),
        	log_mbh,0.5*(merr_log_mbh + perr_log_mbh),mag_sph*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph)
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
        
	logxx = np.arange(-40,40,0.01)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(logxx+np.average(mag_sph),10**yy, color='yellow', ls='-', linewidth=2.)
        m = A[2]
	q = B[2]
	
	AICc_all = akaike.get_AICc_singledataset(mag_sph-np.average(mag_sph),log_mbh,0.5*(merr_log_mbh + perr_log_mbh),m,q,2)
	
       #print 
       #print '-------------------------------'
       #print 'Akaike result coreSersic/Sersic:'
       #print 'Single power law: AICc =', AICc_all
       #print 'Double power law: AICc =', AICc_core
       #print '-------------------------------'
	
       #print 'core-E'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[morph_core=='E_1']-np.average(mag_sph[morph_core=='E_1']),0.5*(perr_mag_sph[morph_core=='E_1']+merr_mag_sph[morph_core=='E_1']),
       #	log_mbh[morph_core=='E_1'],0.5*(merr_log_mbh[morph_core=='E_1'] + perr_log_mbh[morph_core=='E_1']),mag_sph[morph_core=='E_1']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[morph_core=='E_1'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_sph[morph_core=='E_1']),10**yy, color='red', linewidth=2.)
       #m1 = A[2]
       #q1 = B[2]
       #
       #print 'Sersic-bulge'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[BUL_sersic==1]-np.average(mag_sph[BUL_sersic==1]),0.5*(perr_mag_sph[BUL_sersic==1]+merr_mag_sph[BUL_sersic==1]),
       #	log_mbh[BUL_sersic==1],0.5*(merr_log_mbh[BUL_sersic==1] + perr_log_mbh[BUL_sersic==1]),mag_sph[BUL_sersic==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[BUL_sersic==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_sph[BUL_sersic==1]),10**yy, color='blue', ls='--', linewidth=2.)
       #m2 = A[2]
       #q2 = B[2]
       #
       #AICc_ell_bul = akaike.get_AICc_doubledataset(mag_sph[morph_core=='E_1']-np.average(mag_sph[morph_core=='E_1']),log_mbh[morph_core=='E_1'],0.5*(merr_log_mbh[morph_core=='E_1'] + perr_log_mbh[morph_core=='E_1']),m1,q1, 
       #	mag_sph[BUL_sersic==1]-np.average(mag_sph[BUL_sersic==1]),log_mbh[BUL_sersic==1],0.5*(merr_log_mbh[BUL_sersic==1] + perr_log_mbh[BUL_sersic==1]),m2,q2,4)
       #
       #print 
       #print '-------------------------------'
       #print 'Akaike result core ellipticals/sersic bulges:'
       #print 'Double power law: AICc =', AICc_ell_bul
       #print '-------------------------------'
       #
       #cazzata = ELL*[0]
       #for i in range(len(ELL)):
       #	if simplemorphtype[i]=='E' and core[i]==1:
       #		cazzata[i] = 1
       #	elif bulge[i]==1 and core[i]==0:
       #		cazzata[i] = 1
       #		
       #print 'all'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[cazzata==1]-np.average(mag_sph[cazzata==1]),0.5*(perr_mag_sph[cazzata==1]+merr_mag_sph[cazzata==1]),
       #	log_mbh[cazzata==1],0.5*(merr_log_mbh[cazzata==1] + perr_log_mbh[cazzata==1]),mag_sph[cazzata==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[cazzata==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_sph[cazzata==1]),10**yy, color='pink', ls='-', linewidth=2.)
       #m = A[2]
       #q = B[2]
       #
       #AICc_all = akaike.get_AICc_singledataset(mag_sph[cazzata==1]-np.average(mag_sph[cazzata==1]),log_mbh[cazzata==1],0.5*(merr_log_mbh[cazzata==1] + perr_log_mbh[cazzata==1]),m,q,2)
       #
       #print 
       #print '-------------------------------'
       #print 'Akaike result coreSersic ELL and Sersic BUL:'
       #print 'Single power law: AICc =', AICc_all
       #print '-------------------------------'
			
        
       #print 'S0'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']),0.5*(perr_mag_sph[simplemorphtype=='S0']+merr_mag_sph[simplemorphtype=='S0']),
       #	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']),mag_sph[simplemorphtype=='S0']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='S0'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[simplemorphtype=='S0']),10**yy, color='green', linewidth=2.)
       
       #print 'early'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),0.5*(perr_mag_sph[earlytype==1]+merr_mag_sph[earlytype==1]),
       #	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),mag_sph[earlytype==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[earlytype==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-4.,6,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #colorline.colorline(logxx+np.average(mag_sph[earlytype==1]), 10**yy)
       ##ax.plot(logxx+np.average(mag_sph[earlytype==1]),10**yy, color='k', ls='--', linewidth=2.)
       ###### calculates the prediction bands for the given input arrays
       #lpb68_g,upb68_g,logxx_g = predband.predband(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.68,x=logxx[logxx+np.average(mag_sph[earlytype==1])>=-23.5])
       #lpb68_r,upb68_r,logxx_r = predband.predband(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.68,x=logxx[logxx+np.average(mag_sph[earlytype==1])<-23.5])
       ###lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       ##lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       ##### plots a shaded area containing the prediction band
       #ax.fill_between(logxx_g+np.average(mag_sph[earlytype==1]), 10**lpb68_g, 10**upb68_g, alpha=0.15, facecolor='green')
       #ax.fill_between(logxx_r+np.average(mag_sph[earlytype==1]), 10**lpb68_r, 10**upb68_r, alpha=0.15, facecolor='red')
		       
       #print 'S'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='S']-np.average(mag_sph[simplemorphtype=='S']),0.5*(perr_mag_sph[simplemorphtype=='S']+merr_mag_sph[simplemorphtype=='S']),
       #	log_mbh[simplemorphtype=='S'],0.5*(merr_log_mbh[simplemorphtype=='S'] + perr_log_mbh[simplemorphtype=='S']),mag_sph[simplemorphtype=='S']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='S'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_sph[simplemorphtype=='S']),10**yy, color='blue', ls='--', linewidth=4.)
       ###### calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(mag_sph[simplemorphtype=='S']-np.average(mag_sph[simplemorphtype=='S']),log_mbh[simplemorphtype=='S'],A[2],B[2],conf=0.68,x=logxx)
       ###lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       ##lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       ##### plots a shaded area containing the prediction band
       #ax.fill_between(logxx+np.average(mag_sph[simplemorphtype=='S']), 10**lpb68, 10**upb68, alpha=0.15, facecolor='blue')
	
        ax.errorbar(mag_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_mag_sph[simplemorphtype=='merger'],perr_mag_sph[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)

	ax.errorbar(mag_sph[morph_core=='E_1'], mbh[morph_core=='E_1'], 
		xerr=[merr_mag_sph[morph_core=='E_1'],perr_mag_sph[morph_core=='E_1']], 
		yerr=[merr_mbh[morph_core=='E_1'],perr_mbh[morph_core=='E_1']], 
		ecolor='red', marker='o', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mag_sph[morph_core==u'E_0'], mbh[morph_core==u'E_0'], 
		xerr=[merr_mag_sph[morph_core==u'E_0'],perr_mag_sph[morph_core==u'E_0']], 
		yerr=[merr_mbh[morph_core==u'E_0'],perr_mbh[morph_core==u'E_0']], 
		ecolor='red', marker='o', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='S0_1'], mbh[morph_core=='S0_1'], 
		xerr=[merr_mag_sph[morph_core=='S0_1'],perr_mag_sph[morph_core=='S0_1']], 
		yerr=[merr_mbh[morph_core=='S0_1'],perr_mbh[morph_core=='S0_1']], 
		ecolor='green', marker='o', mfc='white', mec='green', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='S0_0'], mbh[morph_core=='S0_0'], 
		xerr=[merr_mag_sph[morph_core=='S0_0'],perr_mag_sph[morph_core=='S0_0']], 
		yerr=[merr_mbh[morph_core=='S0_0'],perr_mbh[morph_core=='S0_0']], 
		ecolor='green', marker='o', mfc='green', mec='green', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='E/S0_1'], mbh[morph_core=='E/S0_1'], 
		xerr=[merr_mag_sph[morph_core=='E/S0_1'],perr_mag_sph[morph_core=='E/S0_1']], 
		yerr=[merr_mbh[morph_core=='E/S0_1'],perr_mbh[morph_core=='E/S0_1']], 
		ecolor='red', marker='*', mfc='white', mec='red', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='E/S0_0'], mbh[morph_core=='E/S0_0'], 
		xerr=[merr_mag_sph[morph_core=='E/S0_0'],perr_mag_sph[morph_core=='E/S0_0']], 
		yerr=[merr_mbh[morph_core=='E/S0_0'],perr_mbh[morph_core=='E/S0_0']], 
		ecolor='red', marker='*', mfc='red', mec='red', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='S_1'], mbh[morph_core=='S_1'], 
        	xerr=[merr_mag_sph[morph_core=='S_1'],perr_mag_sph[morph_core=='S_1']], 
        	yerr=[merr_mbh[morph_core=='S_1'],perr_mbh[morph_core=='S_1']], 
        	ecolor='blue', marker='o', mfc='white', mec='blue', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S_0'], mbh[morph_core=='S_0'], 
        	xerr=[merr_mag_sph[morph_core=='S_0'],perr_mag_sph[morph_core=='S_0']], 
        	yerr=[merr_mbh[morph_core=='S_0'],perr_mbh[morph_core=='S_0']], 
        	ecolor='blue', marker='o', mfc='blue', mec='blue', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S0/S_1'], mbh[morph_core=='S0/S_1'], 
		xerr=[merr_mag_sph[morph_core=='S0/S_1'],perr_mag_sph[morph_core=='S0/S_1']], 
		yerr=[merr_mbh[morph_core=='S0/S_1'],perr_mbh[morph_core=='S0/S_1']], 
		ecolor='blue', marker='*', mfc='white', mec='blue', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S0/S_0'], mbh[morph_core=='S0/S_0'], 
		xerr=[merr_mag_sph[morph_core=='S0/S_0'],perr_mag_sph[morph_core=='S0/S_0']], 
		yerr=[merr_mbh[morph_core=='S0/S_0'],perr_mbh[morph_core=='S0/S_0']], 
		ecolor='blue', marker='*', mfc='blue', mec='blue', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

	ax.set_yscale('log')
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph.pdf', format='pdf', dpi=1000)
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph_morphology.pdf', format='pdf', dpi=1000)
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph_morphology_linregr.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mag_sph_core():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, anc.core_my,  \
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
	gal_id = data[0]
	mbh = data[1].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[2].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[3].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	core = data[4].astype(np.int)
	core_my = data[5].astype(np.int)
	mag_sph = data[6].astype(np.float)
	perr_mag_sph = data[7].astype(np.float)
	merr_mag_sph = data[8].astype(np.float)
	        
        fig, ax = plt.subplots()

        ax.errorbar(mag_sph[core_my==1], mbh[core_my==1], 
		xerr=[merr_mag_sph[core_my==1],perr_mag_sph[core_my==1]], 
		yerr=[merr_mbh[core_my==1],perr_mbh[core_my==1]], 
		ecolor='gray', fmt='wo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[core_my==0], mbh[core_my==0], 
		xerr=[merr_mag_sph[core_my==0],perr_mag_sph[core_my==0]], 
		yerr=[merr_mbh[core_my==0],perr_mbh[core_my==0]], 
		ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
		
	ax.set_yscale('log')
        plt.axis([-15.01,-27.99,10**4.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph_core.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mag_sph_bar():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, anc.bar,  \
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
	gal_id = data[0]
	mbh = data[1].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[2].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[3].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
	core = data[4].astype(np.int)
	bar = data[5].astype(np.int)
	mag_sph = data[6].astype(np.float)
	perr_mag_sph = data[7].astype(np.float)
	merr_mag_sph = data[8].astype(np.float)
	        
        fig, ax = plt.subplots()

        ax.errorbar(mag_sph[bar==1], mbh[bar==1], 
		xerr=[merr_mag_sph[bar==1],perr_mag_sph[bar==1]], 
		yerr=[merr_mbh[bar==1],perr_mbh[bar==1]], 
		ecolor='gray', fmt='w*', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[bar==0], mbh[bar==0], 
		xerr=[merr_mag_sph[bar==0],perr_mag_sph[bar==0]], 
		yerr=[merr_mbh[bar==0],perr_mbh[bar==0]], 
		ecolor='gray', fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
		
	ax.set_yscale('log')
        plt.axis([-18.01,-27.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph_bar.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mag_sph_serind():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, anc.bar,  \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		pysres.log_n_maj_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
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
	bar = data[5].astype(np.int)
	mag_sph = data[6].astype(np.float)
	perr_mag_sph = data[7].astype(np.float)
	merr_mag_sph = data[8].astype(np.float)
	log_n = data[9].astype(np.float)
	n = 10**log_n
	        
        fig, ax = plt.subplots()
	
	#ax.scatter(mag_sph, mbh, c=n)

        ax.errorbar(mag_sph[n>=2], mbh[n>=2], 
        	xerr=[merr_mag_sph[n>=2],perr_mag_sph[n>=2]], 
        	yerr=[merr_mbh[n>=2],perr_mbh[n>=2]], 
        	ecolor='gray', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
      
        ax.errorbar(mag_sph[n<2], mbh[n<2], 
        	xerr=[merr_mag_sph[n<2],perr_mag_sph[n<2]], 
        	yerr=[merr_mbh[n<2],perr_mbh[n<2]], 
        	ecolor='gray', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
		
	ax.set_yscale('log')
        plt.axis([-18.01,-27.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mag_sph_serind.pdf', format='pdf', dpi=1000)
	
def table_mbh_vs_mass_sph_core():

	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.simplemorphtype \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
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
	simplemorphtype = data[8]
	
	mass_sph = 0.5*10**(-0.4*(mag_sph-3.25))
	#perr_mass_sph = 0.5*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) - mass_sph
	#merr_mass_sph = -0.5*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) + mass_sph
	perr_mass_sph = 0.5*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) - mass_sph
	merr_mass_sph = mass_sph - 0.5*10**(-0.4*(mag_sph+perr_mag_sph-3.25))  
	
	log_mass_sph = np.log10(mass_sph)
	perr_log_mass_sph = np.log10(1+perr_mass_sph/mass_sph)
	merr_log_mass_sph = -np.log10(1-merr_mass_sph/mass_sph)
	
	mbh_mass_sph_fileName = 'mbh_mass_sph.dat'
	mbh_mass_sph_file = open(mbh_mass_sph_fileName, 'w')
	
	mbh_mass_sph_file.write('# galaxy  log_mass_sph  +err_log_mass_sph  -err_log_mass_sph  log_mbh  +err_log_mbh  -err_log_mbh  type  core \n')
	
	for i in range(len(mbh)):
		line = gal_id[i] + '   ' + str(log_mass_sph[i]) + '   ' + str(perr_log_mass_sph[i]) + '   ' + str(merr_log_mass_sph[i]) + '   ' 
		line = line + str(log_mbh[i]) + '   ' + str(perr_log_mbh[i]) + '   ' + str(merr_log_mbh[i]) + '   ' 
		line = line + simplemorphtype[i] + '   ' + str(core[i]) + '\n'
		mbh_mass_sph_file.write(line)
	
	mbh_mass_sph_file.close()

	
def mbh_vs_mass_sph_core():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.simplemorphtype, \
		col.color \
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
	simplemorphtype = data[8]
	color = data[9].astype(np.float)
	
	log_ML = 3.98*color+0.13
	ML = 10**log_ML
	
	mass_sph = ML*10**(-0.4*(mag_sph-3.25))
	perr_mass_sph = ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) - mass_sph
	merr_mass_sph = -ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) + mass_sph
	
	log_mass_sph = np.log10(mass_sph)
	perr_log_mass_sph = np.log10(1+perr_mass_sph/mass_sph)
	merr_log_mass_sph = -np.log10(1-merr_mass_sph/mass_sph)
		
	JiangAGNDataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/AGNs-lowmassBHs/GS15-AGN-Jiang.dat'
	JiangAGNDataFile = open(JiangAGNDataFileName)
	
	mass_sph_JiangagnList = []
	mbh_JiangagnList = []
	
	for line in JiangAGNDataFile:
		if line.split()[0] != '#':
			mass_sph_JiangagnList.append(line.split()[3])
			mbh_JiangagnList.append(line.split()[4])
	
	mass_sph_Jiangagn = np.asarray(mass_sph_JiangagnList)
	mbh_Jiangagn = np.asarray(mbh_JiangagnList)	
		
	lowmassAGNDataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/AGNs-lowmassBHs/low-mass.dat'
	lowmassAGNDataFile = open(lowmassAGNDataFileName)
	
	mass_sph_lowmassagnList = []
	mbh_lowmassagnList = []
	
	for line in lowmassAGNDataFile:
		if line.split()[0] != '#':
			mass_sph_lowmassagnList.append(line.split()[2])
			mbh_lowmassagnList.append(line.split()[1])
	
	mass_sph_lowmassagn = np.asarray(mass_sph_lowmassagnList)
	mbh_lowmassagn = np.asarray(mbh_lowmassagnList)	
		
        
        fig, ax = plt.subplots()
	
	ax.errorbar(mass_sph[core==1], mbh[core==1], 
		xerr=[merr_mass_sph[core==1],perr_mass_sph[core==1]], 
		yerr=[merr_mbh[core==1],perr_mbh[core==1]], 
		ecolor='gray', fmt='wo', markersize=10, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mass_sph[core==0], mbh[core==0], 
		xerr=[merr_mass_sph[core==0],perr_mass_sph[core==0]], 
		yerr=[merr_mbh[core==0],perr_mbh[core==0]], 
		ecolor='gray', fmt='ko', markersize=10, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	#ax.scatter(mass_sph_Jiangagn, mbh_Jiangagn, marker='o', color='b', s=15)
	#ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='o', facecolors='none', edgecolors='b', s=80)
	#ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='+', color='b', s=80)
	
        ### fit using FITEXY ###
       #print 'core'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]), 0.5*(perr_log_mass_sph[core==1]+merr_log_mass_sph[core==1]),
       #	log_mbh[core==1], 0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(10**(logxx+np.average(log_mass_sph[core==1])), 10**y_bisec, ls='-.', color='red', linewidth=2.)
       ##label_ell = r'$B_{\rm ell} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-24, 10**6.7, label_ell, fontsize=20)
       

       #print 'core'
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),
       #	0.5*(perr_log_mass_sph[core==1] + merr_log_mass_sph[core==1]),
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),log_mass_sph[core==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[core==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_mass_sph[core==1])),10**yy, color='k', ls='--', linewidth=2.)
       #
       ###### calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.68,x=logxx)
       #lpb95,upb95,logxx = predband.predband(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.95,x=logxx)
       #lpb99,upb99,logxx = predband.predband(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),log_mbh[core==1],A[2],B[2],conf=0.99,x=logxx)
       ##### plots a shaded area containing the prediction band  
       #ax.fill_between(10**(logxx+np.average(log_mass_sph[core==1])), 10**lpb68, 10**upb68, alpha=0.3, facecolor='gray')
       #ax.fill_between(10**(logxx+np.average(log_mass_sph[core==1])), 10**lpb95, 10**upb95, alpha=0.2, facecolor='gray')
       #ax.fill_between(10**(logxx+np.average(log_mass_sph[core==1])), 10**lpb99, 10**upb99, alpha=0.1, facecolor='gray')
       #
       #
       #
       #
       ##print float(mass_sph_Jiangagn)
       ##log_mass_sph_Jiangagn = np.log10(mass_sph_Jiangagn)
       ##log_mass_sph_lowmassagn = np.log10(mass_sph_lowmassagn)
       ##log_mass_sph_sersic = np.concatenate(log_mass_sph[core==0], log_mass_sph_Jiangagn, log_mass_sph_lowmassagn, axis=0)
       ##log_mbh_sersic = np.concatenate((log_mbh[core==0], np.log10(mbh_Jiangagn), np.log10(mbh_lowmassagn)), axis=0)
       #
       #print 'sersic'
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph_sersic-np.average(log_mass_sph_sersic),
       #	(log_mass_sph_sersic*[0.0]+0.1),
       #	log_mbh_sersic,(log_mass_sph_sersic*[0.0]+0.1),log_mass_sph_sersic*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph_sersic)
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_mass_sph_sersic)), 10**yy, color='k', ls='-', linewidth=2.)
	

	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**7.5,10**12.5,10**4,10**11.2])
        plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mass_sph_core.pdf', format='pdf', dpi=1000)
	
	
def mbh_vs_mass_sph():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
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
	
	mass_sph = 0.5*10**(-0.4*(mag_sph-3.25))
	perr_mass_sph = 0.5*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) - mass_sph
	merr_mass_sph = -0.5*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) + mass_sph
	
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
        	if ELL[i]==1 and core[i]==1:
        		ELL_core[i] = 1
        	elif ELL[i]==1 and core[i]==0:
        		ELL_sersic[i] = 1
        	elif ELL[i]==0 and core[i]==1:
        		BUL_core[i] = 1
        	elif ELL[i]==0 and core[i]==0:
        		BUL_sersic[i] = 1
        
        fig, ax = plt.subplots()


        ax.errorbar(mass_sph[simplemorphtype=='E'], mbh[simplemorphtype=='E'], 
		xerr=[merr_mass_sph[simplemorphtype=='E'],perr_mass_sph[simplemorphtype=='E']], 
		yerr=[merr_mbh[simplemorphtype=='E'],perr_mbh[simplemorphtype=='E']], 
		ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'], 
		xerr=[merr_mass_sph[simplemorphtype=='S0'],perr_mass_sph[simplemorphtype=='S0']], 
		yerr=[merr_mbh[simplemorphtype=='S0'],perr_mbh[simplemorphtype=='S0']], 
		ecolor='green', fmt='go', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'], 
		xerr=[merr_mass_sph[simplemorphtype=='E/S0'],perr_mass_sph[simplemorphtype=='E/S0']], 
		yerr=[merr_mbh[simplemorphtype=='E/S0'],perr_mbh[simplemorphtype=='E/S0']], 
		ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[simplemorphtype=='S'], mbh[simplemorphtype=='S'], 
        	xerr=[merr_mass_sph[simplemorphtype=='S'],perr_mass_sph[simplemorphtype=='S']], 
        	yerr=[merr_mbh[simplemorphtype=='S'],perr_mbh[simplemorphtype=='S']], 
        	ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], 
		xerr=[merr_mass_sph[simplemorphtype=='S0/S'],perr_mass_sph[simplemorphtype=='S0/S']], 
		yerr=[merr_mbh[simplemorphtype=='S0/S'],perr_mbh[simplemorphtype=='S0/S']], 
		ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_mass_sph[simplemorphtype=='merger'],perr_mass_sph[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)

	#ax.scatter(mag_sph[gal_id=='m106'], mbh[gal_id=='m106'], color='yellow')
	#ax.scatter(mag_sph[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'], color='g')
	#ax.scatter(mag_sph[simplemorphtype=='S'], mbh[simplemorphtype=='S'], color='b')
	#ax.scatter(mag_sph[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], marker='*', color='b')
	#ax.scatter(mag_sph[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], marker='*', color='b')
	
	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**8.5,10**12.5,10**5.5,10**11.2])
        plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mbh_vs_mass_sph.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mass_tot():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_tot_eq_moffat_comb \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
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
	mag_tot = data[5].astype(np.float)
	
	mass_tot = 0.5*10**(-0.4*(mag_tot-3.25))
	       
        fig, ax = plt.subplots()

	ax.scatter(mass_tot, mbh, color='k')
	
	JiangAGNDataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/AGNs-lowmassBHs/GS15-AGN-Jiang_expanded.dat'
	JiangAGNDataFile = open(JiangAGNDataFileName)
	
	mass_sph_JiangagnList = []
	mbh_JiangagnList = []
	BT_JiangagnList = []
	
	for line in JiangAGNDataFile:
		if line.split()[0] != '#':
			mass_sph_JiangagnList.append(float(line.split()[3]))
			mbh_JiangagnList.append(float(line.split()[4]))
			BT_JiangagnList.append(float(line.split()[5]))
	
	mass_sph_Jiangagn = np.asarray(mass_sph_JiangagnList)
	mbh_Jiangagn = np.asarray(mbh_JiangagnList)	
	BT_Jiangagn = np.asarray(BT_JiangagnList)
	
	mass_tot_Jiangagn = mass_sph_Jiangagn/BT_Jiangagn
	
	ax.scatter(mass_tot_Jiangagn, mbh_Jiangagn, color='b')	
	
	
	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**8.5,10**12.5,10**4.,10**11.2])
        plt.xlabel(r'$M_{\rm *,tot}\rm~[M_\odot]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mass_tot.pdf', format='pdf', dpi=1000)

def mbh_vs_mag_tot():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.simplemorphtype, \
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
	simplemorphtype = data[4]
	mag_sph = data[5].astype(np.float)
	mag_tot = data[6].astype(np.float)
	perr_mag_sph = data[7].astype(np.float)
	merr_mag_sph = data[8].astype(np.float)
	ELL = data[9].astype(np.int)
	
	err_mag_tot = mag_tot*[0.0] + 0.1
	
	earlytype = ELL*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = ELL*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
	
        fig, ax = plt.subplots()

       #print 'core'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_tot[core==1]-np.average(mag_tot[core==1]),err_mag_tot[core==1],
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),mag_tot[core==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_tot[core==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_tot[core==1]),10**yy, color='k', ls='-', linewidth=2.)
       #
       #print 'ser'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_tot[core==0]-np.average(mag_tot[core==0]),err_mag_tot[core==0],
       #	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),mag_tot[core==0]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_tot[core==0])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(logxx+np.average(mag_tot[core==0]),10**yy, color='k', ls='--', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print 'core-Sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_tot[core==1]-np.average(mag_tot[core==1]), err_mag_tot[core==1],
       #	log_mbh[core==1], 0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_tot[core==1]), 10**y_bisec, ls='--', color='b', linewidth=2.)
       ##label_cS = r'$B_{\rm c-Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-19.5, 10**10.5, label_cS, fontsize=20)
       #
       #### fit using FITEXY ###
       #print 'Sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_tot[core==0]-np.average(mag_tot[core==0]), err_mag_tot[core==0],
       #	log_mbh[core==0], 0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]))
       #logxx = np.arange(-40,40,0.01)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(logxx+np.average(mag_tot[core==0]), 10**y_bisec, ls='-', color='b', linewidth=2.)
       ##label_S = r'$B_{\rm Ser} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-19.5, 10**10, label_S, fontsize=20)
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
        
	ax.scatter(mag_tot[simplemorphtype=='E'], mbh[simplemorphtype=='E'], marker='o', color='r', s=40)
        ax.scatter(mag_tot[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'], marker='o', color='g', s=40)
	ax.scatter(mag_tot[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'], marker='*', color='g', s=60)	
	ax.scatter(mag_tot[simplemorphtype=='S'], mbh[simplemorphtype=='S'], marker='o', color='b', s=40)	
	ax.scatter(mag_tot[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], marker='*', color='g', s=60)	
	ax.scatter(mag_tot[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker='d', color='k', s=60)	
	
	
	ax.set_yscale('log')
        plt.axis([-21.01,-27.99,10**5.5,10**11.2])
        plt.xlabel(r'$MAG_{\rm tot}\rm~[mag]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mag_tot.pdf', format='pdf', dpi=1000)
	
def sigma_vs_mag_sph():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.core, anc.simplemorphtype, anc.sigma, \
		pysres.mag_sph_eq_moffat_comb, pysres.mag_tot_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.bar \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	core = data[1].astype(np.int)
	simplemorphtype = data[2]
	sigma = data[3].astype(np.float)
	sigma[gal_id=='n3079'] = 105
	mag_sph = data[4].astype(np.float)
	mag_tot = data[5].astype(np.float)
	perr_mag_sph = data[6].astype(np.float)
	merr_mag_sph = data[7].astype(np.float)
	bar = data[8]
	
	log_sigma = np.log10(sigma)
	err_log_sigma = sigma*[0.0] + np.log10(1.05)
	
	earlytype = core*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
	bulge = core*[0]
	bulge[simplemorphtype=='S0'] = 1
	bulge[simplemorphtype=='S'] = 1
	bulge[simplemorphtype=='S0/S'] = 1
	
        fig, ax = plt.subplots()

       #print 'E'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']),err_mag_sph[simplemorphtype=='E'],
       #	log_sigma[simplemorphtype=='E'],0.5*(perr_log_sigma[simplemorphtype=='E'] + merr_log_sigma[simplemorphtype=='E']),mag_sph[simplemorphtype=='E']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='E'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[simplemorphtype=='E']),yy, color='red', linewidth=2.)
       #
       #print 'S0'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']),err_mag_sph[simplemorphtype=='S0'],
       #	log_sigma[simplemorphtype=='S0'],0.5*(perr_log_sigma[simplemorphtype=='S0'] + merr_log_sigma[simplemorphtype=='S0']),mag_sph[simplemorphtype=='S0']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='S0'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[simplemorphtype=='S0']),yy, color='green', linewidth=2.)
       #
       #print 'early'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),err_mag_sph[earlytype==1],
       #	log_sigma[earlytype==1],0.5*(perr_log_sigma[earlytype==1] + merr_log_sigma[earlytype==1]),mag_sph[earlytype==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[earlytype==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[earlytype==1]),yy, color='k', ls='--', linewidth=2.)
       #
       #print 'bulge'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[bulge==1]-np.average(mag_sph[bulge==1]),err_mag_sph[bulge==1],
       #	log_sigma[bulge==1],0.5*(perr_log_sigma[bulge==1] + merr_log_sigma[bulge==1]),mag_sph[bulge==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[bulge==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[bulge==1]),yy, color='y', ls='-', linewidth=2.)
       #
       #print 'S'
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='S']-np.average(mag_sph[simplemorphtype=='S']),err_mag_sph[simplemorphtype=='S'],
       #	log_sigma[simplemorphtype=='S'],0.5*(perr_log_sigma[simplemorphtype=='S'] + merr_log_sigma[simplemorphtype=='S']),mag_sph[simplemorphtype=='S']*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='S'])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-40,40,0.01)
       #yy = (A[2]*(logxx) + B[2])
       ##ax.plot(logxx+np.average(mag_sph[simplemorphtype=='S']),yy, color='blue', linewidth=2.)
       
	
        ax.scatter(mag_sph[simplemorphtype=='E'], log_sigma[simplemorphtype=='E'], marker='o', color='red', s=40)
        ax.scatter(mag_sph[simplemorphtype=='S0'], log_sigma[simplemorphtype=='S0'], marker='o', color='green', s=40)
        ax.scatter(mag_sph[simplemorphtype=='E/S0'], log_sigma[simplemorphtype=='E/S0'], marker='*', color='green', s=130)
        ax.scatter(mag_sph[simplemorphtype=='S'], log_sigma[simplemorphtype=='S'], marker='o', color='blue', s=40)
        ax.scatter(mag_sph[simplemorphtype=='S0/S'], log_sigma[simplemorphtype=='S0/S'], marker='*', color='green', s=130)
        ax.scatter(mag_sph[simplemorphtype=='merger'], log_sigma[simplemorphtype=='merger'], marker='d', color='k', s=60)
	
	ax.scatter(mag_sph[bar==1], log_sigma[bar==1], marker='*', color='white', s=40)
	
        plt.axis([-19.01,-28.99,1.91,2.59])
        plt.ylabel(r'$\log(\sigma\rm~[km~s^{-1}])$', labelpad=15)
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        plt.show()
	#plt.savefig(path_scalrel_plots + 'mag_sph_vs_sigma_bars.pdf', format='pdf', dpi=1000)
	#plt.savefig(path_scalrel_plots + 'mag_sph_vs_sigma_linregr.pdf', format='pdf', dpi=1000)
	
def mbh_vs_sigma():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.sigma, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH \
		FROM Ancillary AS anc \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	simplemorphtype = data[1]
	sigma = data[2].astype(np.float)
	sigma[gal_id=='n3079'] = 105
	log_sigma = np.log10(sigma)
	err_log_sigma = sigma*[0.0] + np.log10(1.05)
	mbh = data[3].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[4].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_mbh = data[5].astype(np.float)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)
		
	earlytype = sigma*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
	
        fig, ax = plt.subplots()

        print 'all'
        A,B,Aerr,Berr,covAB=bces.bces(log_sigma-np.average(log_sigma),err_log_sigma,
        	log_mbh,0.5*(perr_log_mbh + merr_log_mbh),log_sigma*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_sigma)
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        xx = np.arange(50,500,1)
	logxx = np.log10(xx)
        yy = (A[2]*(logxx-np.average(log_sigma)) + B[2])
        ax.plot(xx,10**yy, color='k', linewidth=2.)
	
        ax.scatter(sigma[simplemorphtype=='E'], mbh[simplemorphtype=='E'], marker='o', color='red', s=40)
        ax.scatter(sigma[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'], marker='o', color='green', s=40)
        ax.scatter(sigma[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'], marker='*', color='green', s=60)
        ax.scatter(sigma[simplemorphtype=='S'], mbh[simplemorphtype=='S'], marker='o', color='blue', s=40)
        ax.scatter(sigma[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], marker='*', color='green', s=60)
        ax.scatter(sigma[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], marker='d', color='k', s=40)
	ax.set_xscale('log')
	ax.set_yscale('log')
	
        plt.axis([80,400,10**5.5,10**10.9])
        plt.xlabel(r'$\sigma\rm~[km~s^{-1}]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH}\rm~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_sigma.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mu_0(axis):
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, \
		pysres.log_n_' + axis + '_moffat_comb, pysres.mu_e_' + axis + '_moffat_comb, \
		errV.perr_mu_0, errV.merr_mu_0, \
		anc.simplemorphtype, anc.core  \
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
	simplemorphtype = data[8]
	core = data[9].astype(np.int)
	
	# compute mu_0
	b = mu_e * [0.0]
	for i in range(len(b)):
		b[i] = b_n.computeb_n(n[i])
	
	mu_0 = mu_e - 2.5*b/np.log(10)

       #A,B,Aerr,Berr,covAB=bces.bces(mu_0-np.average(mu_0),0.5*(merr_mu_0+perr_mu_0),log_mbh,0.5*(merr_log_mbh + perr_log_mbh),mu_0*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mu_0)
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       ##xx = np.arange(0.001,30,0.1)
       #logxx = np.arange(-40,40,0.1)
       ##print A[2], B[2]
       ## calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.68,x=logxx)
       #lpb95,upb95,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       #lpb99,upb99,logxx = predband.predband(mu_0-np.average(mu_0),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       #yy = (A[2]*(logxx) + B[2])
       ##print lpb, upb,xx
	
        fig, ax = plt.subplots()
        ax.errorbar(mu_0[simplemorphtype=='E'], mbh[simplemorphtype=='E'], 
		xerr=[merr_mu_0[simplemorphtype=='E'],perr_mu_0[simplemorphtype=='E']], 
		yerr=[merr_mbh[simplemorphtype=='E'],perr_mbh[simplemorphtype=='E']], 
		ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mu_0[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0'], 
		xerr=[merr_mu_0[simplemorphtype=='E/S0'],perr_mu_0[simplemorphtype=='E/S0']], 
		yerr=[merr_mbh[simplemorphtype=='E/S0'],perr_mbh[simplemorphtype=='E/S0']], 
		ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mu_0[simplemorphtype=='S0'], mbh[simplemorphtype=='S0'], 
		xerr=[merr_mu_0[simplemorphtype=='S0'],perr_mu_0[simplemorphtype=='S0']], 
		yerr=[merr_mbh[simplemorphtype=='S0'],perr_mbh[simplemorphtype=='S0']], 
		ecolor='green', fmt='go', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mu_0[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S'], 
		xerr=[merr_mu_0[simplemorphtype=='S0/S'],perr_mu_0[simplemorphtype=='S0/S']], 
		yerr=[merr_mbh[simplemorphtype=='S0/S'],perr_mbh[simplemorphtype=='S0/S']], 
		ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mu_0[simplemorphtype=='S'], mbh[simplemorphtype=='S'], 
		xerr=[merr_mu_0[simplemorphtype=='S'],perr_mu_0[simplemorphtype=='S']], 
		yerr=[merr_mbh[simplemorphtype=='S'],perr_mbh[simplemorphtype=='S']], 
		ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(mu_0[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_mu_0[simplemorphtype=='merger'],perr_mu_0[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        
	#plt.scatter(logn, mbh, c='black')
        ax.set_yscale('log')
        #xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        #xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        #ax.set_xticks(xticks)
        #ax.set_xticklabels(xticks_labels)
       #ax.plot(logxx+np.average(mu_0),10**yy, color='green', linewidth=2.)
       ## plots a shaded area containing the prediction band  
       #plt.fill_between(logxx+np.average(mu_0), 10**lpb68, 10**upb68, alpha=0.15, facecolor='green')
       #ax.plot(logxx+np.average(mu_0),10**lpb68, color='green')
       #ax.plot(logxx+np.average(mu_0),10**upb68, color='green')
       #plt.fill_between(logxx+np.average(mu_0), 10**lpb95, 10**upb95, alpha=0.1, facecolor='green')
       #ax.plot(logxx+np.average(mu_0),10**lpb95, color='green')
       #ax.plot(logxx+np.average(mu_0),10**upb95, color='green')
       #plt.fill_between(logxx+np.average(mu_0), 10**lpb99, 10**upb99, alpha=0.05, facecolor='green')
       #ax.plot(logxx+np.average(mu_0),10**lpb99, color='green')
       #ax.plot(logxx+np.average(mu_0),10**upb99, color='green')
        plt.axis([-10.1,-34.9,10**5.5,10**11.2])
        plt.xlabel(r'$\mu_0^{\rm ' + axis + r'} \rm ~[mag~arcsec^{-2}]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.17)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mbh_vs_mu_0_' + axis + '.pdf', format='pdf', dpi=1000)
	#plt.clf()
	
		

		
def mu_0_vs_n(axis):
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, \
		pysres.log_n_' + axis + '_moffat_comb, pysres.mu_e_' + axis + '_moffat_comb, \
		errV.perr_log_n, errV.merr_log_n, errV.perr_mu_0, errV.merr_mu_0 \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS pysres ON anc.gal_id = pysres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	simplemorphtype = data[1]
	core = data[2].astype(np.int)
	log_n = data[3].astype(np.float)
	n = 10**log_n
	mu_e = data[4].astype(np.float)
	perr_log_n = data[5].astype(np.float)
	merr_log_n = data[6].astype(np.float)
	perr_mu_0 = data[7].astype(np.float)
	merr_mu_0 = data[8].astype(np.float)
	
	# compute mu_0
	b = mu_e * [0.0]
	for i in range(len(b)):
		b[i] = b_n.computeb_n(n[i])
	
	mu_0 = mu_e - 2.5*b/np.log(10)
	
	earlytype = core*[0]
	earlytype[simplemorphtype=='E'] = 1
	earlytype[simplemorphtype=='S0'] = 1
	earlytype[simplemorphtype=='E/S0'] = 1
		
        fig, ax = plt.subplots()
	
        A,B,Aerr,Berr,covAB=bces.bces(log_n[earlytype==1]-np.average(log_n[earlytype==1]),
        	0.5*(perr_log_n[earlytype==1]+merr_log_n[earlytype==1]),
        	mu_0[earlytype==1],
        	0.5*(merr_mu_0[earlytype==1] + perr_mu_0[earlytype==1]),
        	log_n[earlytype==1]*[0.0])
        print '---------------------------------'
        print 'early'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_n[earlytype==1])
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
        #lpb68,upb68,logxx = predband.predband(log_n-np.average(log_n),mu_0,A[2],B[2],conf=0.68,x=logxx)
        #lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),mu_0,A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),mu_0,A[2],B[2],conf=0.99,x=logxx)
        yy = (A[2]*(logxx) + B[2])
        #print lpb, upb,xx
        #ax.plot(logxx+np.average(log_n[earlytype==1]),yy, color='k', ls='--', linewidth=2.)
       
       #A,B,Aerr,Berr,covAB=bces.bces(log_n[simplemorphtype=='S']-np.average(log_n[simplemorphtype=='S']),
       #	0.5*(perr_log_n[simplemorphtype=='S']+merr_log_n[simplemorphtype=='S']),
       #	log_mbh[simplemorphtype=='S'],
       #	0.5*(merr_log_mbh[simplemorphtype=='S'] + perr_log_mbh[simplemorphtype=='S']),
       #	log_n[simplemorphtype=='S']*[0.0])
       #print '---------------------------------'
       #print 'spi'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_n[simplemorphtype=='S'])
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
       ##lpb68,upb68,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.68,x=logxx)
       ##lpb95,upb95,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.95,x=logxx)
       ##lpb99,upb99,logxx = predband.predband(log_n-np.average(log_n),log_mbh,A[2],B[2],conf=0.99,x=logxx)
       #yy = (A[2]*(logxx) + B[2])
       ##print lpb, upb,xx
       #ax.plot(logxx+np.average(log_n[simplemorphtype=='S']),10**yy, color='b', ls='-', linewidth=2.)
	
        ax.errorbar(log_n[simplemorphtype=='E'], mu_0[simplemorphtype=='E'], xerr=[merr_log_n[simplemorphtype=='E'],perr_log_n[simplemorphtype=='E']], 
        	yerr=[merr_mu_0[simplemorphtype=='E'],perr_mu_0[simplemorphtype=='E']], ecolor='red', fmt='ro', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[simplemorphtype=='E/S0'], mu_0[simplemorphtype=='E/S0'], xerr=[merr_log_n[simplemorphtype=='E/S0'],perr_log_n[simplemorphtype=='E/S0']], 
        	yerr=[merr_mu_0[simplemorphtype=='E/S0'],perr_mu_0[simplemorphtype=='E/S0']], ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[simplemorphtype=='S0'], mu_0[simplemorphtype=='S0'], xerr=[merr_log_n[simplemorphtype=='S0'],perr_log_n[simplemorphtype=='S0']], 
        	yerr=[merr_mu_0[simplemorphtype=='S0'],perr_mu_0[simplemorphtype=='S0']], ecolor='green', fmt='go', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[simplemorphtype=='S'], mu_0[simplemorphtype=='S'], xerr=[merr_log_n[simplemorphtype=='S'],perr_log_n[simplemorphtype=='S']], 
        	yerr=[merr_mu_0[simplemorphtype=='S'],perr_mu_0[simplemorphtype=='S']], ecolor='blue', fmt='bo', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[simplemorphtype=='S0/S'], mu_0[simplemorphtype=='S0/S'], xerr=[merr_log_n[simplemorphtype=='S0/S'],perr_log_n[simplemorphtype=='S0/S']], 
        	yerr=[merr_mu_0[simplemorphtype=='S0/S'],perr_mu_0[simplemorphtype=='S0/S']], ecolor='green', fmt='g*', markersize=20, elinewidth=1.2, capthick=1.2, barsabove=False) 
        ax.errorbar(log_n[simplemorphtype=='merger'], mu_0[simplemorphtype=='merger'], xerr=[merr_log_n[simplemorphtype=='merger'],perr_log_n[simplemorphtype=='merger']], 
        	yerr=[merr_mu_0[simplemorphtype=='merger'],perr_mu_0[simplemorphtype=='merger']], ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 

        xticks = np.log10(np.asarray([0.5,1,2,3,4,5,6,10]))
        xticks_labels = ['$0.5$','$1$','$2$','$3$','$4$','$5$','$6$','$10$']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels)
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
        plt.axis([np.log10(0.25),np.log10(18),-10,-35])
        plt.xlabel(r'$n^{\rm ' + axis + '}$', labelpad=15)
        plt.ylabel(r'$\mu_0^{\rm ' + axis + r'} \rm ~[mag~arcsec^{-2}]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15)
        #plt.show()
	plt.savefig(path_scalrel_plots + 'mu_0_vs_n_' + axis + '.pdf', format='pdf', dpi=1000)
	


def main():
	#mbh_vs_n('maj')	
	#mbh_vs_n('eq')		
	#mbh_vs_mu_0('maj')
	#mbh_vs_mu_0('eq')
	#mbh_vs_mag_sph()
	#mbh_vs_mag_sph_core()
	#mbh_vs_mag_sph_bar()
	#mbh_vs_mag_sph_serind()
	#mbh_vs_mag_tot()
	#mbh_vs_mass_tot()
	#mbh_vs_logn('eq')		
	#mbh_vs_logr_e('comparison')
	#mbh_vs_logr_e()
	#mag_sph_vs_n('maj')
	#mag_sph_vs_n('eq')
	#KMAG_sph_vs_mag_sph()
	#log_mbh_vs_KMAG_sph()
	#sani_vs_me()
	#BTratios_GS13_vs_me()
	#mag_tot_GS13_vs_me()
	#sigma_vs_mag_sph()
	#mbh_vs_sigma()
	#mu_0_vs_n('maj')
	#mu_0_vs_n('eq')
	#mbh_vs_mass_sph()
	mbh_vs_mass_sph_core()
	#table_mbh_vs_mass_sph_core()
	
main()
