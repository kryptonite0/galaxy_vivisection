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
from instruments.markers import markers

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 22})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/MbhMsph/images/'


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


def mag_lit_vs_mag_my():
	
	connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id,  \
		physres.mag_sph_eq_moffat_comb, \
		physsgs.mag_sph, \
		physsani.mag_sph, \
		physlasker.mag_sph, \
		physvika.mag_sph \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN KmagSGS2013PhysicalUnits AS physsgs ON anc.gal_id = physsgs.gal_id \
		JOIN LiteratureDecompositionsSanietal2011PhysicalUnits AS physsani ON anc.gal_id = physsani.gal_id \
		JOIN LiteratureDecompositionsLaskeretal2014PhysicalUnits AS physlasker ON anc.gal_id = physlasker.gal_id \
		JOIN LiteratureDecompositionsVikaetal2012PhysicalUnits AS physvika ON anc.gal_id = physvika.gal_id \
		WHERE anc.fit1D_done = 1;'
	cur.execute(getdata_query)
        datalist = cur.fetchall()
        data= np.asarray(datalist).transpose()
	gal_id = data[0]
	mag_sph_my = data[1].astype(np.float)	
	mag_sph_sgs = data[2].astype(np.float)	
	mag_sph_sani = data[3].astype(np.float)	
	mag_sph_lasker = data[4].astype(np.float)	
	mag_sph_vika = data[5].astype(np.float)	

	fig, ax = plt.subplots()
	ax.scatter(mag_sph_my, mag_sph_sgs+0.27, c='r', s=60)
	ax.scatter(mag_sph_my, mag_sph_sani, c='g', s=60, marker='s')
	ax.scatter(mag_sph_my, mag_sph_lasker+0.27, c='orange', s=120, marker='*')
	ax.scatter(mag_sph_my, mag_sph_vika+0.27, c='yellow', s=80, marker='^')
	ax.plot(np.arange(-28,-18,1),np.arange(-28,-18,1), ls='--', color='k', lw=2)
	ax.axis([-19.01,-27.99,-19.01,-27.99])
	ax.set_xlabel(r'$MAG_{\rm sph} \rm ~[mag]~(this~work)$', labelpad=15)
	ax.set_ylabel(r'$MAG_{\rm sph} \rm ~[mag]~(literature)$', labelpad=15)
	
	## legend
	ax.scatter([-19.6],[-27.2], c='r', s=60)
	ax.scatter([-19.6],[-26.5], c='g', s=60, marker='s')
	ax.scatter([-19.6],[-25.8], c='orange', s=120, marker='*')
	ax.scatter([-19.6],[-25.1], c='yellow', s=80, marker='^')
	
	ax.text(-19.9, -27., 'SGS13', fontsize=20)
	ax.text(-19.9, -26.3, 'S+11', fontsize=20)
	ax.text(-19.9, -25.6, 'L+14', fontsize=20)
	ax.text(-19.9, -24.9, 'V+12', fontsize=20)
	
        plt.subplots_adjust(left=0.15,bottom=0.15)
	#plt.show()
	plt.savefig(path_paper_figures + 'mag_lit_vs_mag_my.pdf', format='pdf', dpi=1000)
	
def mbh_vs_mass_sph_agn():
	
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
	perr_mass_sph = ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) - mass_sph
	merr_mass_sph = -ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) + mass_sph
	
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
	
	JiangAGNDataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/AGNs-lowmassBHs/GS15-AGN-Jiang_expanded.dat'
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

        print 'BCES all'
        print 'n', len(log_mass_sph)
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph-np.average(log_mass_sph),
        	0.5*(perr_log_mass_sph + merr_log_mass_sph),
        	log_mbh,0.5*(merr_log_mbh + perr_log_mbh),log_mass_sph*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_mass_sph-np.average(log_mass_sph), log_mbh, B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_mass_sph-np.average(log_mass_sph), log_mbh, B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_mass_sph-np.average(log_mass_sph), log_mbh, B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_mass_sph-np.average(log_mass_sph), log_mbh, B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph)
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(10**(logxx+np.average(log_mass_sph)),10**yy, color='k', ls='-', linewidth=2.)
        #colorline.colorline(10**(logxx+np.average(log_mass_sph)), 10**yy, cmap=green_red)
       
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
        			
        ax.fill_between(10**(logxx+np.average(log_mass_sph)), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='k')
       
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY all'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph-np.average(log_mass_sph),
       #	0.5*(perr_log_mass_sph + merr_log_mass_sph),
       #	log_mbh,0.5*(merr_log_mbh + perr_log_mbh))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       ############################
        
        print 'BCES early'
        print 'n', len(log_mass_sph[earlytype==1])
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),
        	0.5*(perr_log_mass_sph[earlytype==1] + merr_log_mass_sph[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),log_mass_sph[earlytype==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]), log_mbh[earlytype==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]), log_mbh[earlytype==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]), log_mbh[earlytype==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]), log_mbh[earlytype==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph[earlytype==1])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
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
       
        ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='r')
       
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY early'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),
       #	0.5*(perr_log_mass_sph[earlytype==1] + merr_log_mass_sph[earlytype==1]),
       #	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #############################
       #
        print 'sersic bul of spi'
        print 'n', len(log_mass_sph[morph_core=='S_0'])
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),
        	0.5*(perr_log_mass_sph[morph_core=='S_0'] + merr_log_mass_sph[morph_core=='S_0']),
        	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']),log_mass_sph[morph_core=='S_0']*[0.0])
        absscat_0 = absolutescatter.get_absscatter(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph[morph_core=='S_0'])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(10**(logxx+np.average(log_mass_sph[morph_core=='S_0'])),10**yy, color='b', ls='-', linewidth=2.)
       
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
        			
        ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='S_0'])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='b')
       
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY sersic bul of spi'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),
       #	0.5*(perr_log_mass_sph[morph_core=='S_0'] + merr_log_mass_sph[morph_core=='S_0']),
       #	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #print 'ellipticals'
       #print 'n', len(log_mass_sph[simplemorphtype=='E'])
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']),
       #	0.5*(perr_log_mass_sph[simplemorphtype=='E'] + merr_log_mass_sph[simplemorphtype=='E']),
       #	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E']),log_mass_sph[simplemorphtype=='E']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[simplemorphtype=='E'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY ellipticals'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[simplemorphtype=='E']-np.average(log_mass_sph[simplemorphtype=='E']),
       #	0.5*(perr_log_mass_sph[simplemorphtype=='E'] + merr_log_mass_sph[simplemorphtype=='E']),
       #	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print 'S0s'
       #print 'n', len(log_mass_sph[simplemorphtype=='S0'])
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']),
       #	0.5*(perr_log_mass_sph[simplemorphtype=='S0'] + merr_log_mass_sph[simplemorphtype=='S0']),
       #	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']),log_mass_sph[simplemorphtype=='S0']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[simplemorphtype=='S0'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY S0s'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[simplemorphtype=='S0']-np.average(log_mass_sph[simplemorphtype=='S0']),
       #	0.5*(perr_log_mass_sph[simplemorphtype=='S0'] + merr_log_mass_sph[simplemorphtype=='S0']),
       #	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print 'core'
       #print 'n', len(log_mass_sph[core==1])
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),
       #	0.5*(perr_log_mass_sph[core==1] + merr_log_mass_sph[core==1]),
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),log_mass_sph[core==1]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]), log_mbh[core==1], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]), log_mbh[core==1], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]), log_mbh[core==1], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]), log_mbh[core==1], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[core==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(log_mass_sph[core==1])),10**yy, color='k', ls='-.', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY core'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[core==1]-np.average(log_mass_sph[core==1]),
       #	0.5*(perr_log_mass_sph[core==1] + merr_log_mass_sph[core==1]),
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #print 'sersic'
       #print 'n', len(log_mass_sph[core==0])
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]),
       #	0.5*(perr_log_mass_sph[core==0] + merr_log_mass_sph[core==0]),
       #	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),log_mass_sph[core==0]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]), log_mbh[core==0], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]), log_mbh[core==0], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]), log_mbh[core==0], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]), log_mbh[core==0], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[core==0])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(log_mass_sph[core==0])),10**yy, color='k', ls=':', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[core==0]-np.average(log_mass_sph[core==0]),
       #	0.5*(perr_log_mass_sph[core==0] + merr_log_mass_sph[core==0]),
       #	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]))
       #print '-----------------------------'
       #print '-----------------------------'

	###################################

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

	ax.scatter(mass_sph_Jiangagn, mbh_Jiangagn, marker='o', color='b', s=15)
	ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='o', facecolors='none', edgecolors='b', s=80)
	ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='+', color='b', s=80)

	ax.errorbar(mass_sph[morph_core=='E_1'], mbh[morph_core=='E_1'], 
		xerr=[merr_mass_sph[morph_core=='E_1'],perr_mass_sph[morph_core=='E_1']], 
		yerr=[merr_mbh[morph_core=='E_1'],perr_mbh[morph_core=='E_1']], 
		ecolor='red', marker='o', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mass_sph[morph_core==u'E_0'], mbh[morph_core==u'E_0'], 
		xerr=[merr_mass_sph[morph_core==u'E_0'],perr_mass_sph[morph_core==u'E_0']], 
		yerr=[merr_mbh[morph_core==u'E_0'],perr_mbh[morph_core==u'E_0']], 
		ecolor='red', marker='o', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[morph_core=='S0_1'], mbh[morph_core=='S0_1'], 
		xerr=[merr_mass_sph[morph_core=='S0_1'],perr_mass_sph[morph_core=='S0_1']], 
		yerr=[merr_mbh[morph_core=='S0_1'],perr_mbh[morph_core=='S0_1']], 
		ecolor='red', marker='^', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[morph_core=='S0_0'], mbh[morph_core=='S0_0'], 
		xerr=[merr_mass_sph[morph_core=='S0_0'],perr_mass_sph[morph_core=='S0_0']], 
		yerr=[merr_mbh[morph_core=='S0_0'],perr_mbh[morph_core=='S0_0']], 
		ecolor='red', marker='^', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[morph_core=='E/S0_1'], mbh[morph_core=='E/S0_1'], 
		xerr=[merr_mass_sph[morph_core=='E/S0_1'],perr_mass_sph[morph_core=='E/S0_1']], 
		yerr=[merr_mbh[morph_core=='E/S0_1'],perr_mbh[morph_core=='E/S0_1']], 
		ecolor='red', marker='*', mfc='white', mec='red', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[morph_core=='E/S0_0'], mbh[morph_core=='E/S0_0'], 
		xerr=[merr_mass_sph[morph_core=='E/S0_0'],perr_mass_sph[morph_core=='E/S0_0']], 
		yerr=[merr_mbh[morph_core=='E/S0_0'],perr_mbh[morph_core=='E/S0_0']], 
		ecolor='red', marker='*', mfc='red', mec='red', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mass_sph[morph_core=='S_1'], mbh[morph_core=='S_1'], 
        	xerr=[merr_mass_sph[morph_core=='S_1'],perr_mass_sph[morph_core=='S_1']], 
        	yerr=[merr_mbh[morph_core=='S_1'],perr_mbh[morph_core=='S_1']], 
        	ecolor='blue', marker='s', mfc='white', mec='blue', mew=1.5, markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[morph_core=='S_0'], mbh[morph_core=='S_0'], 
        	xerr=[merr_mass_sph[morph_core=='S_0'],perr_mass_sph[morph_core=='S_0']], 
        	yerr=[merr_mbh[morph_core=='S_0'],perr_mbh[morph_core=='S_0']], 
        	ecolor='blue', marker='s', mfc='blue', mec='blue', markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[morph_core=='S0/S_1'], mbh[morph_core=='S0/S_1'], 
		xerr=[merr_mass_sph[morph_core=='S0/S_1'],perr_mass_sph[morph_core=='S0/S_1']], 
		yerr=[merr_mbh[morph_core=='S0/S_1'],perr_mbh[morph_core=='S0/S_1']], 
		ecolor='blue', marker='v', mfc='white', mec='blue', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[morph_core=='S0/S_0'], mbh[morph_core=='S0/S_0'], 
		xerr=[merr_mass_sph[morph_core=='S0/S_0'],perr_mass_sph[morph_core=='S0/S_0']], 
		yerr=[merr_mbh[morph_core=='S0/S_0'],perr_mbh[morph_core=='S0/S_0']], 
		ecolor='blue', marker='v', mfc='blue', mec='blue', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mass_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_mass_sph[simplemorphtype=='merger'],perr_mass_sph[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mass_sph[morph_core=='outlier'], mbh[morph_core=='outlier'], 
		xerr=[merr_mass_sph[morph_core=='outlier'],perr_mass_sph[morph_core=='outlier']], 
		yerr=[merr_mbh[morph_core=='outlier'],perr_mbh[morph_core=='outlier']], 
		ecolor='gray', fmt='kx', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.scatter(mass_sph[morph_core=='outlier'], mbh[morph_core=='outlier'], marker='x', c='k', s=100, lw=2, **scatter_kwargs)
		
	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**7.9,10**12.3,10**4.1,10**11.2])
        plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mass_sph_agn.pdf', format='pdf', dpi=1000)

def mbh_vs_mag_sph_psb():

	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		anc.ELLIPTICAL_my, anc.simplemorphtype, \
		pysres.log_n_eq_moffat_comb, \
		anc.pseudobulge_KH13 \
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
	cbar.ax.set_ylabel(r'$n_{\rm sph}$', rotation=90)
	
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
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
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
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
	print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
        print '---------------------------------'
        
	slope_n_grt_2 = A[2]
	intercept_n_grt_2 = B[2]
	
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot((logxx+np.average(mag_sph[n>2])),10**yy, color='k', ls='-', linewidth=2.)
	
        ### fit using FITEXY ###
        print '-----------------------------'
        print 'FITEXY n>2'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[n>2]-np.average(mag_sph[n>2]),
        	0.5*(perr_mag_sph[n>2] + merr_mag_sph[n>2]),
        	log_mbh[n>2],0.5*(merr_log_mbh[n>2] + perr_log_mbh[n>2]))
        print '-----------------------------'
        print '-----------------------------'
	
	###########################

	# make inset
	ins = plt.axes([.62, .25, .18, .2])
	ins.axis([0.01,11,-1.99,1.99])
	ins.xaxis.set_ticks([1,5,10])
	ins.yaxis.set_ticks([-1,0,1])
	ins.set_xlabel(r'$n_{\rm sph}$')
	ins.set_ylabel(r'offset')
	ins.plot([-1,15], [0,0], ls='-', c='k')
	ins.plot([2,2], [-3,3], ls='--', c='k')
	offset = log_mbh - (slope_n_grt_2*(mag_sph-np.average(mag_sph[n>2])) + intercept_n_grt_2)
	ins.scatter(n, offset, marker='o', facecolor='k', edgecolor='k', s=15, **scatter_kwargs)
	#ins.scatter(n[pseudobulge_KH13==1], offset[pseudobulge_KH13==1], marker='o', facecolor='b', edgecolor='gray', s=20, **scatter_kwargs)
	
	ax.set_yscale('log')
        ax.axis([-19.01,-27.99,10**5.1,10**11.2])
        ax.set_xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        ax.set_ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.99,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mag_sph_psb.pdf', format='pdf', dpi=1000)

def mbh_vs_mag_tot():

	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	lowerlimit = [u'm94', u'n3079', u'n4388', u'n4945']
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.core, \
		physres.mag_tot_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		res.delta_eq_moffat_comb, \
		anc.ELLIPTICAL_my, anc.simplemorphtype, \
		col.color  \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN OneDFitResults AS res ON anc.gal_id = res.gal_id \
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
	mag_tot = data[5].astype(np.float)
	perr_mag_sph = data[6].astype(np.float)
	merr_mag_sph = data[7].astype(np.float)
	delta = data[8].astype(np.float)
	ELL = data[9].astype(np.int)
	simplemorphtype = data[10]
	color = data[11].astype(np.float)
	
	err_mag_tot = mag_tot*[0.0] + 0.2
	       
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
        	if gal_id[i] not in lowerlimit:
			morph_coreList.append(simplemorphtype[i] + '_' + str(core[i]))
		else:
			morph_coreList.append('lowerlimit')	
        morph_core = np.asarray(morph_coreList)
	
	all = ELL*[0]
	for i in range(len(ELL)):
        	if gal_id[i] not in lowerlimit:
			all[i] = 1

        fig, ax = plt.subplots()
	
	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

       #	
       #print 'BCES all'
       #print 'n', len(mag_tot[all==1])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_tot[all==1]-np.average(mag_tot[all==1]),
       #	err_mag_tot[all==1],
       #	log_mbh[all==1],0.5*(merr_log_mbh[all==1] + perr_log_mbh[all==1]),mag_tot[all==1]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_tot[all==1]-np.average(mag_tot[all==1]), log_mbh[all==1], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_tot[all==1]-np.average(mag_tot[all==1]), log_mbh[all==1], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_tot[all==1]-np.average(mag_tot[all==1]), log_mbh[all==1], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_tot[all==1]-np.average(mag_tot[all==1]), log_mbh[all==1], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_tot[all==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_tot[all==1])),10**yy, color='k', ls='-', linewidth=2.)
       ##colorline.colorline(10**(logxx+np.average(mag_tot[all==1])), 10**yy, cmap=green_red)
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
       #ax.fill_between((logxx+np.average(mag_tot[all==1])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='k')
       #
       #### fit using FITEXY ###
       #print 'FITEXY all'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_tot[all==1]-np.average(mag_tot[all==1]),
       #       err_mag_tot[all==1],
       #       log_mbh[all==1],0.5*(merr_log_mbh[all==1] + perr_log_mbh[all==1]))
       #
       #print '----------------------------'
       #print '----------------------------'
       #
       #######################################
       #
       #print 'BCES early'
       #print 'n', len(mag_tot[earlytype==1])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]),
       #	err_mag_tot[earlytype==1],
       #	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),mag_tot[earlytype==1]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]), log_mbh[earlytype==1], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]), log_mbh[earlytype==1], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]), log_mbh[earlytype==1], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]), log_mbh[earlytype==1], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_tot[earlytype==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_tot[earlytype==1])),10**yy, color='r', ls='--', linewidth=2.)
       ##colorline.colorline(10**(logxx+np.average(mag_tot[earlytype==1])), 10**yy, cmap=green_red)
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
       #ax.fill_between((logxx+np.average(mag_tot[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='r')
	
       #### fit using FITEXY ###
       #print 'FITEXY early'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_tot[earlytype==1]-np.average(mag_tot[earlytype==1]),
       #       err_mag_tot[earlytype==1],
       #       log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]))
       #
       #print '----------------------------'
       #print '----------------------------'
       #
       #######################################
	
        print 'sersic spirals'
        print 'n', len(mag_tot[morph_core=='S_0'])
        print stats.spearmanr(mag_tot[morph_core=='S_0'], log_mbh[morph_core=='S_0'])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']),
       #	err_mag_tot[morph_core=='S_0'],
       #	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']),mag_tot[morph_core=='S_0']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_tot[morph_core=='S_0'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-50,20,1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_tot[morph_core=='S_0'])),10**yy, color='r', ls='--', linewidth=2.)
       ##colorline.colorline(10**(logxx+np.average(mag_tot[morph_core=='S_0'])), 10**yy, cmap=green_red)
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
       #ax.fill_between((logxx+np.average(mag_tot[morph_core=='S_0'])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='b')
       #
       #### fit using FITEXY ###
       #print 'FITEXY sersic spirals'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_tot[morph_core=='S_0']-np.average(mag_tot[morph_core=='S_0']),
       #	err_mag_tot[morph_core=='S_0'],
       #	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']))
       #
       #print '----------------------------'
       #print '----------------------------'
	
	######################################
	
	
	ax.errorbar(mag_tot[morph_core=='E_1'], mbh[morph_core=='E_1'], 
		xerr=[err_mag_tot[morph_core=='E_1'],err_mag_tot[morph_core=='E_1']], 
		yerr=[merr_mbh[morph_core=='E_1'],perr_mbh[morph_core=='E_1']], 
		ecolor='red', marker='o', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mag_tot[morph_core==u'E_0'], mbh[morph_core==u'E_0'], 
		xerr=[err_mag_tot[morph_core==u'E_0'],err_mag_tot[morph_core==u'E_0']], 
		yerr=[merr_mbh[morph_core==u'E_0'],perr_mbh[morph_core==u'E_0']], 
		ecolor='red', marker='o', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_tot[morph_core=='S0_1'], mbh[morph_core=='S0_1'], 
		xerr=[err_mag_tot[morph_core=='S0_1'],err_mag_tot[morph_core=='S0_1']], 
		yerr=[merr_mbh[morph_core=='S0_1'],perr_mbh[morph_core=='S0_1']], 
		ecolor='red', marker='^', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_tot[morph_core=='S0_0'], mbh[morph_core=='S0_0'], 
		xerr=[err_mag_tot[morph_core=='S0_0'],err_mag_tot[morph_core=='S0_0']], 
		yerr=[merr_mbh[morph_core=='S0_0'],perr_mbh[morph_core=='S0_0']], 
		ecolor='red', marker='^', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_tot[morph_core=='E/S0_1'], mbh[morph_core=='E/S0_1'], 
		xerr=[err_mag_tot[morph_core=='E/S0_1'],err_mag_tot[morph_core=='E/S0_1']], 
		yerr=[merr_mbh[morph_core=='E/S0_1'],perr_mbh[morph_core=='E/S0_1']], 
		ecolor='red', marker='*', mfc='white', mec='red', mew=1.5, markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_tot[morph_core=='E/S0_0'], mbh[morph_core=='E/S0_0'], 
		xerr=[err_mag_tot[morph_core=='E/S0_0'],err_mag_tot[morph_core=='E/S0_0']], 
		yerr=[merr_mbh[morph_core=='E/S0_0'],perr_mbh[morph_core=='E/S0_0']], 
		ecolor='red', marker='*', mfc='red', mec='red', markersize=20, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_tot[morph_core=='S_1'], mbh[morph_core=='S_1'], 
        	xerr=[err_mag_tot[morph_core=='S_1'],err_mag_tot[morph_core=='S_1']], 
        	yerr=[merr_mbh[morph_core=='S_1'],perr_mbh[morph_core=='S_1']], 
        	ecolor='blue', marker='s', mfc='white', mec='blue', mew=1.5, markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_tot[morph_core=='S_0'], mbh[morph_core=='S_0'], 
        	xerr=[err_mag_tot[morph_core=='S_0'],err_mag_tot[morph_core=='S_0']], 
        	yerr=[merr_mbh[morph_core=='S_0'],perr_mbh[morph_core=='S_0']], 
        	ecolor='blue', marker='s', mfc='blue', mec='blue', markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_tot[morph_core=='S0/S_1'], mbh[morph_core=='S0/S_1'], 
		xerr=[err_mag_tot[morph_core=='S0/S_1'],err_mag_tot[morph_core=='S0/S_1']], 
		yerr=[merr_mbh[morph_core=='S0/S_1'],perr_mbh[morph_core=='S0/S_1']], 
		ecolor='blue', marker='v', mfc='white', mec='blue', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_tot[morph_core=='S0/S_0'], mbh[morph_core=='S0/S_0'], 
		xerr=[err_mag_tot[morph_core=='S0/S_0'],err_mag_tot[morph_core=='S0/S_0']], 
		yerr=[merr_mbh[morph_core=='S0/S_0'],perr_mbh[morph_core=='S0/S_0']], 
		ecolor='blue', marker='v', mfc='blue', mec='blue', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_tot[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[err_mag_tot[simplemorphtype=='merger'],err_mag_tot[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
 	ax.errorbar(mag_tot[morph_core=='lowerlimit'], mbh[morph_core=='lowerlimit'], 
 	       xerr=-0.1, 
 	       yerr=[merr_mbh[morph_core=='lowerlimit'],perr_mbh[morph_core=='lowerlimit']], 
 	       xuplims=True, ecolor='blue', elinewidth=1.2, capthick=1.2, ls=' ', barsabove=False)
	
	#ax.scatter(mag_tot[morph_core=='lowerlimit'], mbh[morph_core=='lowerlimit'], marker='>', c='b', s=100, lw=2, **scatter_kwargs)
		
	ax.set_yscale('log')
        plt.axis([-21.5,-27.99,10**5.1,10**11.2])
        plt.xlabel(r'$MAG_{\rm gal}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mag_tot.pdf', format='pdf', dpi=1000)

def mbh_vs_mag_sph():
	
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

        print 'BCES all'
        print 'n', len(mag_sph)
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph-np.average(mag_sph),
        	0.5*(perr_mag_sph + merr_mag_sph),
        	log_mbh,0.5*(merr_log_mbh + perr_log_mbh),mag_sph*[0.0])
        absscat_0 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(mag_sph-np.average(mag_sph), log_mbh, B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph)
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot((logxx+np.average(mag_sph)),10**yy, color='k', ls='-', linewidth=2.)
        #colorline.colorline(10**(logxx+np.average(mag_sph)), 10**yy, cmap=green_red)
       
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
        			
        ax.fill_between((logxx+np.average(mag_sph)), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='k')
       
        ### fit using FITEXY ###
        print '-----------------------------'
        print 'FITEXY all'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph-np.average(mag_sph),
        	0.5*(perr_mag_sph + merr_mag_sph),
        	log_mbh,0.5*(merr_log_mbh + perr_log_mbh))
        print '-----------------------------'
        print '-----------------------------'
       
        ###########################
       
        print 'BCES early'
        print 'n', len(mag_sph[earlytype==1])
        A,B,Aerr,Berr,covAB=bces.bces(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),
        	0.5*(perr_mag_sph[earlytype==1] + merr_mag_sph[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]),mag_sph[earlytype==1]*[0.0])
        absscat_0 = absolutescatter.get_absscatter(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]), log_mbh[earlytype==1], B[0], A[0])
        absscat_1 = absolutescatter.get_absscatter(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]), log_mbh[earlytype==1], B[1], A[1])
        absscat_2 = absolutescatter.get_absscatter(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]), log_mbh[earlytype==1], B[2], A[2])
        absscat_3 = absolutescatter.get_absscatter(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]), log_mbh[earlytype==1], B[3], A[3])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(mag_sph[earlytype==1])
        print
        print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
        print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot((logxx+np.average(mag_sph[earlytype==1])),10**yy, color='r', ls='--', linewidth=2.)
        #colorline.colorline(10**(logxx+np.average(mag_sph[earlytype==1])), 10**yy, cmap=green_red)
       
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
        			
        ax.fill_between((logxx+np.average(mag_sph[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='r')
       
        ### fit using FITEXY ###
        print '-----------------------------'
        print 'FITEXY early'
        A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[earlytype==1]-np.average(mag_sph[earlytype==1]),
        	0.5*(perr_mag_sph[earlytype==1] + merr_mag_sph[earlytype==1]),
        	log_mbh[earlytype==1],0.5*(merr_log_mbh[earlytype==1] + perr_log_mbh[earlytype==1]))
        print '-----------------------------'
        print '-----------------------------'
       
       #############################
       #
       #print 'sersic bul of spi'
       #print 'n', len(mag_sph[morph_core=='S_0'])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']),
       #	0.5*(perr_mag_sph[morph_core=='S_0'] + merr_mag_sph[morph_core=='S_0']),
       #	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']),mag_sph[morph_core=='S_0']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']), log_mbh[morph_core=='S_0'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[morph_core=='S_0'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_sph[morph_core=='S_0'])),10**yy, color='b', ls='-', linewidth=2.)
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
       #ax.fill_between((logxx+np.average(mag_sph[morph_core=='S_0'])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='b')
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY sersic bul of spi'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[morph_core=='S_0']-np.average(mag_sph[morph_core=='S_0']),
       #	0.5*(perr_mag_sph[morph_core=='S_0'] + merr_mag_sph[morph_core=='S_0']),
       #	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #print 'ellipticals'
       #print 'n', len(mag_sph[simplemorphtype=='E'])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']),
       #	0.5*(perr_mag_sph[simplemorphtype=='E'] + merr_mag_sph[simplemorphtype=='E']),
       #	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E']),mag_sph[simplemorphtype=='E']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']), log_mbh[simplemorphtype=='E'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='E'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY ellipticals'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[simplemorphtype=='E']-np.average(mag_sph[simplemorphtype=='E']),
       #	0.5*(perr_mag_sph[simplemorphtype=='E'] + merr_mag_sph[simplemorphtype=='E']),
       #	log_mbh[simplemorphtype=='E'],0.5*(merr_log_mbh[simplemorphtype=='E'] + perr_log_mbh[simplemorphtype=='E']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print 'S0s'
       #print 'n', len(mag_sph[simplemorphtype=='S0'])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']),
       #	0.5*(perr_mag_sph[simplemorphtype=='S0'] + merr_mag_sph[simplemorphtype=='S0']),
       #	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']),mag_sph[simplemorphtype=='S0']*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']), log_mbh[simplemorphtype=='S0'], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[simplemorphtype=='S0'])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY S0s'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[simplemorphtype=='S0']-np.average(mag_sph[simplemorphtype=='S0']),
       #	0.5*(perr_mag_sph[simplemorphtype=='S0'] + merr_mag_sph[simplemorphtype=='S0']),
       #	log_mbh[simplemorphtype=='S0'],0.5*(merr_log_mbh[simplemorphtype=='S0'] + perr_log_mbh[simplemorphtype=='S0']))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #
       ###########################
       #
       #print 'core'
       #print 'n', len(mag_sph[core==1])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==1]-np.average(mag_sph[core==1]),
       #	0.5*(perr_mag_sph[core==1] + merr_mag_sph[core==1]),
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]),mag_sph[core==1]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_sph[core==1]-np.average(mag_sph[core==1]), log_mbh[core==1], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_sph[core==1]-np.average(mag_sph[core==1]), log_mbh[core==1], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph[core==1]-np.average(mag_sph[core==1]), log_mbh[core==1], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_sph[core==1]-np.average(mag_sph[core==1]), log_mbh[core==1], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[core==1])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_sph[core==1])),10**yy, color='k', ls='-.', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY core'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==1]-np.average(mag_sph[core==1]),
       #	0.5*(perr_mag_sph[core==1] + merr_mag_sph[core==1]),
       #	log_mbh[core==1],0.5*(merr_log_mbh[core==1] + perr_log_mbh[core==1]))
       #print '-----------------------------'
       #print '-----------------------------'
       #
       #print 'sersic'
       #print 'n', len(mag_sph[core==0])
       #A,B,Aerr,Berr,covAB=bces.bces(mag_sph[core==0]-np.average(mag_sph[core==0]),
       #	0.5*(perr_mag_sph[core==0] + merr_mag_sph[core==0]),
       #	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]),mag_sph[core==0]*[0.0])
       #absscat_0 = absolutescatter.get_absscatter(mag_sph[core==0]-np.average(mag_sph[core==0]), log_mbh[core==0], B[0], A[0])
       #absscat_1 = absolutescatter.get_absscatter(mag_sph[core==0]-np.average(mag_sph[core==0]), log_mbh[core==0], B[1], A[1])
       #absscat_2 = absolutescatter.get_absscatter(mag_sph[core==0]-np.average(mag_sph[core==0]), log_mbh[core==0], B[2], A[2])
       #absscat_3 = absolutescatter.get_absscatter(mag_sph[core==0]-np.average(mag_sph[core==0]), log_mbh[core==0], B[3], A[3])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(mag_sph[core==0])
       #print
       #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0]), 'Delta =', "{0:.4f}".format(absscat_0)
       #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1]), 'Delta =', "{0:.4f}".format(absscat_1)
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2]), 'Delta =', "{0:.4f}".format(absscat_2)
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3]), 'Delta =', "{0:.4f}".format(absscat_3)
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot((logxx+np.average(mag_sph[core==0])),10**yy, color='k', ls=':', linewidth=2.)
       #
       #### fit using FITEXY ###
       #print '-----------------------------'
       #print 'FITEXY sersic'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(mag_sph[core==0]-np.average(mag_sph[core==0]),
       #	0.5*(perr_mag_sph[core==0] + merr_mag_sph[core==0]),
       #	log_mbh[core==0],0.5*(merr_log_mbh[core==0] + perr_log_mbh[core==0]))
       #print '-----------------------------'
       #print '-----------------------------'

	###################################
	
		
	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

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
		ecolor='red', marker='^', mfc='white', mec='red', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
        ax.errorbar(mag_sph[morph_core=='S0_0'], mbh[morph_core=='S0_0'], 
		xerr=[merr_mag_sph[morph_core=='S0_0'],perr_mag_sph[morph_core=='S0_0']], 
		yerr=[merr_mbh[morph_core=='S0_0'],perr_mbh[morph_core=='S0_0']], 
		ecolor='red', marker='^', mfc='red', mec='red', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)
	
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
        	ecolor='blue', marker='s', mfc='white', mec='blue', mew=1.5, markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S_0'], mbh[morph_core=='S_0'], 
        	xerr=[merr_mag_sph[morph_core=='S_0'],perr_mag_sph[morph_core=='S_0']], 
        	yerr=[merr_mbh[morph_core=='S_0'],perr_mbh[morph_core=='S_0']], 
        	ecolor='blue', marker='s', mfc='blue', mec='blue', markersize=9, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S0/S_1'], mbh[morph_core=='S0/S_1'], 
		xerr=[merr_mag_sph[morph_core=='S0/S_1'],perr_mag_sph[morph_core=='S0/S_1']], 
		yerr=[merr_mbh[morph_core=='S0/S_1'],perr_mbh[morph_core=='S0/S_1']], 
		ecolor='blue', marker='v', mfc='white', mec='blue', mew=1.5, markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[morph_core=='S0/S_0'], mbh[morph_core=='S0/S_0'], 
		xerr=[merr_mag_sph[morph_core=='S0/S_0'],perr_mag_sph[morph_core=='S0/S_0']], 
		yerr=[merr_mbh[morph_core=='S0/S_0'],perr_mbh[morph_core=='S0/S_0']], 
		ecolor='blue', marker='v', mfc='blue', mec='blue', markersize=12, ls=' ', elinewidth=1.2, capthick=1.2, barsabove=False)

        ax.errorbar(mag_sph[simplemorphtype=='merger'], mbh[simplemorphtype=='merger'], 
		xerr=[merr_mag_sph[simplemorphtype=='merger'],perr_mag_sph[simplemorphtype=='merger']], 
		yerr=[merr_mbh[simplemorphtype=='merger'],perr_mbh[simplemorphtype=='merger']], 
		ecolor='k', fmt='kd', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.errorbar(mag_sph[morph_core=='outlier'], mbh[morph_core=='outlier'], 
		xerr=[merr_mag_sph[morph_core=='outlier'],perr_mag_sph[morph_core=='outlier']], 
		yerr=[merr_mbh[morph_core=='outlier'],perr_mbh[morph_core=='outlier']], 
		ecolor='gray', fmt='kx', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False)
	
	ax.scatter(mag_sph[morph_core=='outlier'], mbh[morph_core=='outlier'], marker='x', c='k', s=100, lw=2, **scatter_kwargs)
	
	
	
	ax.set_yscale('log')
        plt.axis([-19.01,-27.99,10**5.1,10**11.2])
        plt.xlabel(r'$MAG_{\rm sph}\rm~[mag]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'mbh_vs_mag_sph.pdf', format='pdf', dpi=1000)

def mbh_vs_mass_sph_galsymb_agn():
	
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
	perr_mass_sph = ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) - mass_sph
	merr_mass_sph = -ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) + mass_sph
	
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
	
	JiangAGNDataFileName = '/Users/gsavorgnan/galaxy_vivisection/data/Alister-data/AGNs-lowmassBHs/GS15-AGN-Jiang_expanded.dat'
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

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

	ax.scatter(mass_sph_Jiangagn, mbh_Jiangagn, marker='o', color='b', s=15, alpha=0.5)
	ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='o', facecolors='none', edgecolors='b', s=80, alpha=0.5)
	ax.scatter(mass_sph_lowmassagn, mbh_lowmassagn, marker='+', color='b', s=80, alpha=0.5)

	
	for x0,y0 in zip(mass_sph[simplemorphtype=='E'], mbh[simplemorphtype=='E']):
		markers.elliptical(ax, 'red', np.log10(x0), np.log10(y0), 0.08)
	
 	for x0,y0 in zip(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
		markers.lenticular(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.08)
       
	for x0,y0 in zip(mass_sph[simplemorphtype=='S0'], mbh[simplemorphtype=='S0']):
		markers.lenticular(ax, 'darkorange', np.log10(x0), np.log10(y0), 0.08)
	
        for x0,y0 in zip(mass_sph[simplemorphtype=='S0/S'], mbh[simplemorphtype=='S0/S']):
		markers.spiral(ax, 'blue', np.log10(x0), np.log10(y0), 0.04)
			
        for x0,y0 in zip(mass_sph[simplemorphtype=='S'], mbh[simplemorphtype=='S']):
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
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
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
				
	ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='r')
	
	##### calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.68,x=logxx)
        #lpb95,upb95,logxx = predband.predband(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_mass_sph[earlytype==1]-np.average(log_mass_sph[earlytype==1]),log_mbh[earlytype==1],A[2],B[2],conf=0.99,x=logxx)
        #### plots a shaded area containing the prediction band  
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**lpb68, 10**upb68, alpha=0.1, facecolor='r')
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**lpb95, 10**upb95, alpha=0.07, facecolor='r')
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[earlytype==1])), 10**lpb99, 10**upb99, alpha=0.04, facecolor='r')
	
       #print 'sersic bul'
       #print 'n', len(log_mass_sph[BUL_sersic==1])
       #A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[BUL_sersic==1]-np.average(log_mass_sph[BUL_sersic==1]),
       #	0.5*(perr_log_mass_sph[BUL_sersic==1] + merr_log_mass_sph[BUL_sersic==1]),
       #	log_mbh[BUL_sersic==1],0.5*(merr_log_mbh[BUL_sersic==1] + perr_log_mbh[BUL_sersic==1]),log_mass_sph[BUL_sersic==1]*[0.0])
       #print '---------------------------------'
       #print 'y = A*(x-<x>) + B '
       #print '<x> =', np.average(log_mass_sph[BUL_sersic==1])
       #print
       ##print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
       ##print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
       #print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
       ##print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
       #print '---------------------------------'
       #
       #logxx = np.arange(-10,20,0.1)
       #yy = (A[2]*(logxx) + B[2])
       #ax.plot(10**(logxx+np.average(log_mass_sph[BUL_sersic==1])),10**yy, color='b', ls='-', linewidth=2.)
       #
       ###### calculates the prediction bands for the given input arrays
       #lpb68,upb68,logxx = predband.predband(log_mass_sph[BUL_sersic==1]-np.average(log_mass_sph[BUL_sersic==1]),log_mbh[BUL_sersic==1],A[2],B[2],conf=0.68,x=logxx)
       #lpb95,upb95,logxx = predband.predband(log_mass_sph[BUL_sersic==1]-np.average(log_mass_sph[BUL_sersic==1]),log_mbh[BUL_sersic==1],A[2],B[2],conf=0.95,x=logxx)
       #lpb99,upb99,logxx = predband.predband(log_mass_sph[BUL_sersic==1]-np.average(log_mass_sph[BUL_sersic==1]),log_mbh[BUL_sersic==1],A[2],B[2],conf=0.99,x=logxx)
       ##### plots a shaded area containing the prediction band  
       #ax.fill_between(10**(logxx+np.average(log_mass_sph[BUL_sersic==1])), 10**lpb68, 10**upb68, alpha=0.1, facecolor='b')
       ##ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='E_1'])), 10**lpb95, 10**upb95, alpha=0.07, facecolor='r')
       ##ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='E_1'])), 10**lpb99, 10**upb99, alpha=0.04, facecolor='r')
        
        print 'sersic bul of spi'
        print 'n', len(log_mass_sph[morph_core=='S_0'])
        A,B,Aerr,Berr,covAB=bces.bces(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),
		0.5*(perr_log_mass_sph[morph_core=='S_0'] + merr_log_mass_sph[morph_core=='S_0']),
        	log_mbh[morph_core=='S_0'],0.5*(merr_log_mbh[morph_core=='S_0'] + perr_log_mbh[morph_core=='S_0']),log_mass_sph[morph_core=='S_0']*[0.0])
        print '---------------------------------'
        print 'y = A*(x-<x>) + B '
        print '<x> =', np.average(log_mass_sph[morph_core=='S_0'])
        print
        #print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
        #print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
        print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
        #print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
        print '---------------------------------'
       
        logxx = np.arange(-10,20,0.1)
        yy = (A[2]*(logxx) + B[2])
        ax.plot(10**(logxx+np.average(log_mass_sph[morph_core=='S_0'])),10**yy, color='b', ls='-', linewidth=2.)
	
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
				
	ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='S_0'])), 10**yy_lo, 10**yy_up, alpha=0.1, facecolor='b')
	
	#print gal_id[morph_core=='S_0'],log_mbh[morph_core=='S_0']
	##### calculates the prediction bands for the given input arrays
        #lpb68,upb68,logxx = predband.predband(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),log_mbh[morph_core=='S_0'],A[2],B[2],conf=0.68,x=logxx)
        #lpb95,upb95,logxx = predband.predband(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),log_mbh[morph_core=='S_0'],A[2],B[2],conf=0.95,x=logxx)
        #lpb99,upb99,logxx = predband.predband(log_mass_sph[morph_core=='S_0']-np.average(log_mass_sph[morph_core=='S_0']),log_mbh[morph_core=='S_0'],A[2],B[2],conf=0.99,x=logxx)
        #### plots a shaded area containing the prediction band  
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='S_0'])), 10**lpb68, 10**upb68, alpha=0.1, facecolor='b')
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='E_1'])), 10**lpb95, 10**upb95, alpha=0.07, facecolor='r')
        #ax.fill_between(10**(logxx+np.average(log_mass_sph[morph_core=='E_1'])), 10**lpb99, 10**upb99, alpha=0.04, facecolor='r')
	
       #### fit using FITEXY ###
       #print 'sersic bul'
       #A,perr_A,merr_A,B,perr_B,merr_B = fitexy.bisect_modfitexy(log_mass_sph[BUL_sersic==1]-np.average(log_mass_sph[BUL_sersic==1]),
       #	0.5*(perr_log_mass_sph[BUL_sersic==1] + merr_log_mass_sph[BUL_sersic==1]),
       #	log_mbh[BUL_sersic==1],0.5*(merr_log_mbh[BUL_sersic==1] + perr_log_mbh[BUL_sersic==1]))
       #logxx = np.arange(-10,20,0.1)
       ## plot bisector relation
       #y_bisec = A + B*logxx
       #ax.plot(10**(logxx+np.average(log_mass_sph[BUL_sersic==1])), 10**y_bisec, ls='-.', color='blue', linewidth=2.)
       ##label_ell = r'$B_{\rm ell} = ' + str("{0:.2f}".format(B)) + r'^{+' + str("{0:.2f}".format(perr_B)) + r'}_{-' + str("{0:.2f}".format(merr_B)) +r'}$'
       ##ax.text(-24, 10**6.7, label_ell, fontsize=20)
	
	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**7.9,10**12.3,10**4.1,10**11.2])
        plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=15)
        plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=15)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        #plt.show()
	plt.savefig(path_paper_figures + 'mbh_vs_mass_sph_galsymb_agn.pdf', format='pdf', dpi=1000)

	
		
		
def main():
	#mag_lit_vs_mag_my()
	mbh_vs_mass_sph_agn()
	#mbh_vs_mag_sph_psb()
	#mbh_vs_mag_tot()
	#mbh_vs_mag_sph()
	#mbh_vs_mass_sph_galsymb_agn()

main()		
		
	
