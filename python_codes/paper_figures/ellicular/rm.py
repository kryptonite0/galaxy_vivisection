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

path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/ellicular/images/'

def Re_vs_mass_sph():
	
	#outliers = [u'n1374', u'n3842exp', u'n4889']
	outliers = []
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, \
		pysres.mag_sph_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		pysres.log_r_e_maj_moffat_comb, \
		errV.perr_log_r_e, errV.merr_log_r_e, \
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
	mag_sph = data[1].astype(np.float)
	perr_mag_sph = data[2].astype(np.float)
	merr_mag_sph = data[3].astype(np.float)
	log_r_e_sph = data[4].astype(np.float)
	perr_log_r_e_sph = data[5].astype(np.float)
	merr_log_r_e_sph = data[6].astype(np.float)
	color = data[7].astype(np.float)
	
	log_ML = 3.98*color+0.13 # meidt+2014
	ML = 10**log_ML
	
	mass_sph = ML*10**(-0.4*(mag_sph-3.25))
	merr_mass_sph = -ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) + mass_sph
	perr_mass_sph = +ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) - mass_sph
	
	r_e_sph = 10**log_r_e_sph
	
	log_mass_sph = np.log10(mass_sph)
	perr_log_mass_sph = np.log10(1+perr_mass_sph/mass_sph)
	merr_log_mass_sph = -np.log10(1-merr_mass_sph/mass_sph)
	
	ok = mag_sph*[0.0]
	ok[gal_id=='n1332'] = 1
	ok[gal_id=='n3115'] = 1
	ok[gal_id=='n3377'] = 1
	ok[gal_id=='n0821'] = 1	
	ok[gal_id=='n4697'] = 1	
	
	mass_sph_n1277 = 2.69*10**11
	r_e_sph_n1277 = 2.12 # kpc
	
	mass_sph_n1271 = 9.17*10**10
	r_e_sph_n1271 = 1.3 # kpc
	
	mass_sph_m1216 = 2.27*10**11
	r_e_sph_m1216 = 3.01 # kpc
	       
        ##################################################################
	
	fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}
	
	ax.scatter(mass_sph[ok==1], r_e_sph[ok==1], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	ax.scatter([mass_sph_n1277], [r_e_sph_n1277], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	ax.scatter([mass_sph_n1271], [r_e_sph_n1271], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	ax.scatter([mass_sph_m1216], [r_e_sph_m1216], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
	

       #for x0,y0 in zip(mass_sph[earlytype==1], mbh[earlytype==1]):
       #	markers.elliptical(ax, 'k', np.log10(x0), np.log10(y0), 0.05)
       #for x0,y0 in zip(mass_sph[gal_id=='n3115'], mbh[gal_id=='n3115']):
       #	markers.elliptical(ax, 'green', np.log10(x0), np.log10(y0), 0.05)	
       #for x0,y0 in zip(mass_sph[gal_id=='n4697'], mbh[gal_id=='n4697']):
       #	markers.elliptical(ax, 'green', np.log10(x0), np.log10(y0), 0.05)	
       #for x0,y0 in zip(mass_sph[gal_id=='n0821'], mbh[gal_id=='n0821']):
       #	markers.elliptical(ax, 'green', np.log10(x0), np.log10(y0), 0.05)	
       #for x0,y0 in zip(mass_sph[gal_id=='n3377'], mbh[gal_id=='n3377']):
       #	markers.elliptical(ax, 'green', np.log10(x0), np.log10(y0), 0.05)	
       #
       ##for x0,y0 in zip(mass_sph[simplemorphtype=='E/S0'], mbh[simplemorphtype=='E/S0']):
       ##	markers.lenticular(ax, 'red', np.log10(x0), np.log10(y0), 0.08)
       #ax.scatter(mass_sph[gal_id=='n1332'], mbh[gal_id=='n1332'], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
       #
       #mass_sph_n1277 = 2.69*10**11
       #mass_sph_n1277_b = mass_sph_n1277/11.65*6	
       #mass_sph_n1277_old = 2.88*10**10
       #mbh_n1277 = 1.7*10**10
       #ax.scatter([mass_sph_n1277], [mbh_n1277], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
       ##ax.scatter([mass_sph_n1277_b], [mbh_n1277], marker=r'$\star$', s=100, color='red', **scatter_kwargs)
       #ax.scatter([mass_sph_n1277_old], [mbh_n1277], marker=r'$\star$', s=500, color='gray', **scatter_kwargs)
       #ax.plot([mass_sph_n1277_old,mass_sph_n1277], [mbh_n1277,mbh_n1277], color='gray', lw=2, ls='--')
       #
       #mass_sph_n1271 = 9.17*10**10
       #mass_sph_n1271_old = 5.4*10**10
       #mbh_n1271 = 3*10**9
       #ax.scatter([mass_sph_n1271], [mbh_n1271], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
       ##ax.plot([mass_sph_n1271_old,mass_sph_n1271], [mbh_n1271,mbh_n1271], color='blue', lw=3)
       #
       #mass_sph_m1216 = 2.27*10**11
       #mbh_m1216 = 10**10 ## upper limit
       ##ax.scatter([mass_sph_m1216], [mbh_m1216], marker='.', s=500, color='red', **scatter_kwargs)
       #ax.scatter([mass_sph_m1216], [mbh_m1216/1.5], marker=r'$\downarrow$', s=500, color='red', **scatter_kwargs)
       #
       #mass_sph_n1332 = mass_sph[gal_id=='n1332']
       #mass_sph_n1332_old = mass_sph_n1332/0.95*0.43
       #mbh_n1332 = mbh[gal_id=='n1332']
       ##ax.plot([mass_sph_n1332_old,mass_sph_n1332], [mbh_n1332,mbh_n1332], color='blue', lw=3)
       #
       #
       #x0 = mass_sph[gal_id=='n3115']
       #y0 = mbh[gal_id=='n3115']
       #ax.text(x0/1.5, 1.4*y0, 'N3115', size=12, color='green')
       #
       #x0 = mass_sph[gal_id=='n4697']
       #y0 = mbh[gal_id=='n4697']
       #ax.text(x0/1.4, y0/2., 'N4697', size=12, color='green')
       #
       #x0 = mass_sph[gal_id=='n0821']
       #y0 = mbh[gal_id=='n0821']
       #ax.text(1.2*x0, y0/1.3, 'N0821', size=12, color='green')
       #
       #x0 = mass_sph[gal_id=='n3377']
       #y0 = mbh[gal_id=='n3377']
       #ax.text(1.2*x0, y0/1.4, 'N3377', size=12, color='green')
       #
       #x0 = mass_sph[gal_id=='n1332']
       #y0 = mbh[gal_id=='n1332']
       #ax.text(x0*1.2, y0/1.3, 'N1332', size=12, color='red')
       #
       #x0 = mass_sph_n1277
       #y0 = mbh_n1277
       #ax.text(x0*1.2, y0/1.3, 'N1277', size=12, color='red')
       #
       #x0 = mass_sph_n1271
       #y0 = mbh_n1271
       #ax.text(x0*1.1, y0*1.3, 'N1271', size=12, color='red')
       #
       #x0 = mass_sph_m1216
       #y0 = mbh_m1216
       #ax.text(x0/1.9, 0.8*y0, 'M1216', size=12, color='red')
       #
       #x0 = float(mass_sph[gal_id=='n4291'])
       #y0 = float(mbh[gal_id=='n4291'])
       #ax.scatter([x0], [y0], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
       #ax.text(x0/2, y0/1.6, 'N4291', size=12, color='red')
       #
       ##x0 = float(mass_sph[gal_id=='n3998'])
       ##y0 = float(mbh[gal_id=='n3998'])
       ##ax.text(x0/0.95, 1.3*y0, 'N3998', size=12, color='k')
	
			
	
	# legend
       #markers.elliptical(ax, 'k', 8.85, 10.9, 0.05)
       #ax.text(10**9.1, 10**10.75, 'Elliptical')
       #ax.scatter([10**8.85], [10**10.45], marker=r'$\star$', s=500, color='red', **scatter_kwargs)
       ##markers.lenticular(ax, 'red', 8.85, 10.45, 0.08)
       #ax.text(10**9.1, 10**10.3, 'Ellicular')
       #markers.lenticular(ax, 'k', 8.85, 10., 0.05)
       #ax.text(10**9.1, 10**9.85, 'Lenticular')
       ##markers.spiral(ax, 'darkorange', 9.8, 10.9, 0.04)
       ##ax.text(10**10.05, 10**10.75, 'S0/Sp')
       ##markers.spiral(ax, 'blue', 9.8, 10.45, 0.04)
       ##ax.text(10**10.05, 10**10.3, 'Sp')
       ##ax.scatter([10**9.8], [10**10.], marker=r'$\star$', s=500, color='k', **scatter_kwargs)	
       ##ax.text(10**10.05, 10**9.85, 'merger')
	
 	ax.set_xscale('log')
	ax.set_yscale('log')
        plt.axis([10**9.8,10**12.1,0.5,50])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
	plt.xlabel(r'Spheroid stellar mass $\rm [M_\odot]$', labelpad=13)
	plt.ylabel(r'Spheroid size $\rm [kpc]$', labelpad=13)
	plt.subplots_adjust(left=0.15,bottom=0.15,right=0.97,top=0.9)
        plt.show()
	#plt.savefig(path_paper_figures + 'rm.pdf', format='pdf', dpi=1000)



Re_vs_mass_sph()