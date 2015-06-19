import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#import matplotlib.colors as mcolors
#import matplotlib.image as image

#from scipy import stats
#from instruments.linear_regression import bces
#from instruments.linear_regression import predband
#from instruments.linear_regression import fitexy
from instruments import b_n
#from instruments.linear_regression import colorline
#from instruments.linear_regression import akaike
#from instruments.linear_regression import absolutescatter
#from instruments.linear_regression import lnr
from instruments.markers import markers_lin, markers

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 26})

linregr1varAllFileName = '/Users/gsavorgnan/galaxy_vivisection/python_codes/BHFP/BHFP_all_1varfit.out'
linregr2varEarlyFileName = '/Users/gsavorgnan/galaxy_vivisection/python_codes/BHFP/BHFP_early_2varfit.out'

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

#path_figures = 

connection = sql3.connect(dbname)
cur = connection.cursor()

getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.bar, \
	anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.sigma, \
	physres.mag_sph_eq_moffat_comb, physres.mag_tot_eq_moffat_comb, \
	physres.log_n_maj_moffat_comb, physres.log_n_eq_moffat_comb, \
	physres.mu_e_maj_moffat_comb, physres.mu_e_eq_moffat_comb, \
	physres.log_r_e_maj_moffat_comb, physres.log_r_e_eq_moffat_comb, \
	errV.perr_mag_sph, errV.merr_mag_sph, \
	errV.perr_log_n, errV.merr_log_n, \
	errV.perr_mu_e, errV.merr_mu_e, \
        errV.perr_log_r_e, errV.merr_log_r_e, \
	errV.perr_mu_0, errV.merr_mu_0 \
	FROM Ancillary AS anc \
	JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
	JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
	WHERE anc.fit1D_done = 1;'

cur.execute(getdata_query)
datalist = cur.fetchall()
data= np.asarray(datalist).transpose()

gal_id = data[0]
simplemorphtype = data[1]
core = data[2].astype(np.int)
bar = data[3].astype(np.int)

earlytype = np.zeros(len(gal_id))
for i in range(len(gal_id)):
    if simplemorphtype[i]=='E' or simplemorphtype[i]=='E/S0' or simplemorphtype[i]=='S0':
        earlytype[i] = 1
latetype = np.zeros(len(gal_id))
for i in range(len(gal_id)):
    if simplemorphtype[i]=='Sp':
        latetype[i] = 1

mbh = data[4].astype(np.float)
log_mbh = np.log10(mbh)
perr_mbh = data[5].astype(np.float)
merr_mbh = data[6].astype(np.float)
perr_log_mbh = np.log10(1 + perr_mbh/mbh)
merr_log_mbh = -np.log10(1 - merr_mbh/mbh)

sigma = data[7].astype(np.float)
# assign value to n3079
sigma[gal_id=='n3079'] = 105
log_sigma = np.log10(sigma)
err_log_sigma = sigma*[0.0] + np.log10(1.05)

mag_sph = data[8].astype(np.float)
mag_tot = data[9].astype(np.float)
perr_mag_sph = data[16].astype(np.float)
merr_mag_sph = data[17].astype(np.float)
err_mag_tot = np.zeros(len(gal_id)) + 0.25

log_n_maj = data[10].astype(np.float)
n_maj = 10**log_n_maj
log_n_eq = data[11].astype(np.float)
n_eq = 10**log_n_eq
perr_log_n = data[18].astype(np.float)
merr_log_n = data[19].astype(np.float)

mu_e_maj = data[12].astype(np.float)
mu_e_eq = data[13].astype(np.float)
perr_mu_e = data[20].astype(np.float)
merr_mu_e = data[21].astype(np.float)

# compute mu_0
b_maj = np.zeros(len(gal_id))
for i in range(len(b_maj)):
	b_maj[i] = b_n.computeb_n(n_maj[i])
b_eq = np.zeros(len(gal_id))
for i in range(len(b_eq)):
	b_eq[i] = b_n.computeb_n(n_eq[i])

mu_0_maj = mu_e_maj - 2.5*b_maj/np.log(10)
mu_0_eq = mu_e_eq - 2.5*b_eq/np.log(10)
perr_mu_0 = data[24].astype(np.float)
merr_mu_0 = data[25].astype(np.float)

log_r_e_maj = data[14].astype(np.float)
#r_e_maj = 10**log_r_e_maj
log_r_e_eq = data[15].astype(np.float)
#r_e_eq = 10**log_r_e_eq
perr_log_r_e = data[22].astype(np.float)
merr_log_r_e = data[23].astype(np.float)

print mu_e_eq


def MAGsph_sigma_Req():
	linregr2varEarlyFile = open(linregr2varEarlyFileName)
	linregr2varEarly = linregr2varEarlyFile.readlines()
	i = linregr2varEarly.index('mag_sph - log_sigma - log_r_e_eq\n')
	log_sigma_mean = float((linregr2varEarly[i+4]).split()[2])
	log_r_e_eq_mean = float((linregr2varEarly[i+5]).split()[2])
	mag_sph_mean = float((linregr2varEarly[i+6]).split()[2])
	A = float((linregr2varEarly[i+7]).split()[2])
	B0 = float((linregr2varEarly[i+8]).split()[2])
	B1 = float((linregr2varEarly[i+9]).split()[2])
	
	X = A+ B0*(log_sigma-log_sigma_mean) + B1*(log_r_e_eq-log_r_e_eq_mean) + mag_sph_mean
	
        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

        for x0,y0 in zip(X[simplemorphtype=='E'], mag_sph[simplemorphtype=='E']):
        	markers_lin.elliptical(ax, 'red', (x0), (y0), 0.08, 0.08)
       
        for x0,y0 in zip(X[simplemorphtype=='E/S0'], mag_sph[simplemorphtype=='E/S0']):
        	markers_lin.lenticular(ax, 'red', x0, y0, 0.14, 0.1)
       
        for x0,y0 in zip(X[simplemorphtype=='S0'], mag_sph[simplemorphtype=='S0']):
        	markers_lin.lenticular(ax, 'darkorange', x0, y0, 0.14, 0.1)
       
        for x0,y0 in zip(X[simplemorphtype=='S0/Sp'], mag_sph[simplemorphtype=='S0/Sp']):
        	markers_lin.spiral(ax, 'darkorange', x0, y0, 0.08, 0.08)
        		
        for x0,y0 in zip(X[simplemorphtype=='Sp'], mag_sph[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.08, 0.08)
        	
        for x0,y0 in zip(X[simplemorphtype=='Sp'], mag_sph[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.08, 0.08)
       
        ax.scatter(X[simplemorphtype=='merger'], mag_sph[simplemorphtype=='merger'], marker=r'$\star$', s=300, color='k', **scatter_kwargs)	
 
        plt.xlabel(r'$\alpha \log(\sigma~{\rm [km~s^{-1}]}) + \beta \log(R_{\rm e}~{\rm [kpc]})$', labelpad=13)
        plt.ylabel(r'$MAG_{\rm sph} \rm ~[mag]$', labelpad=13)
	plt.subplots_adjust(left=0.17,bottom=0.17,right=0.97,top=0.9)
        #plt.axis([-13.5, -5.01, -29, -19])
	
	#ax.text(-13, -20, r'\alpha = ' + )
	
        plt.show()
	#plt.savefig('MAGsph_sigma_Re.pdf', format='pdf', dpi=1000)

	linregr2varEarlyFile.close()
	
def Req_sigma_Ie():
	linregr2varEarlyFile = open(linregr2varEarlyFileName)
	linregr2varEarly = linregr2varEarlyFile.readlines()
	i = linregr2varEarly.index('log_r_e_eq - log_sigma - mu_e_eq\n')
	log_sigma_mean = float((linregr2varEarly[i+4]).split()[2])
	mu_e_eq_mean = float((linregr2varEarly[i+5]).split()[2])
	log_r_e_eq_mean = float((linregr2varEarly[i+6]).split()[2])
	A = float((linregr2varEarly[i+7]).split()[2])
	B0 = float((linregr2varEarly[i+8]).split()[2])
	B1 = float((linregr2varEarly[i+9]).split()[2])
	
	X = B0*(log_sigma) + (B1)*(mu_e_eq) 
	print B1/B0
	
        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

       #for x0,y0 in zip(X[simplemorphtype=='E'], log_r_e_eq[simplemorphtype=='E']):
       #	markers_lin.elliptical(ax, 'red', (x0), (y0), 0.08, 0.08)
       #
       #for x0,y0 in zip(X[simplemorphtype=='E/S0'], log_r_e_eq[simplemorphtype=='E/S0']):
       #	markers_lin.lenticular(ax, 'red', x0, y0, 0.14, 0.1)
       #
       #for x0,y0 in zip(X[simplemorphtype=='S0'], log_r_e_eq[simplemorphtype=='S0']):
       #	markers_lin.lenticular(ax, 'darkorange', x0, y0, 0.14, 0.1)
       #
       #for x0,y0 in zip(X[simplemorphtype=='S0/Sp'], log_r_e_eq[simplemorphtype=='S0/Sp']):
       #	markers_lin.spiral(ax, 'darkorange', x0, y0, 0.08, 0.08)
       #		
       #for x0,y0 in zip(X[simplemorphtype=='Sp'], log_r_e_eq[simplemorphtype=='Sp']):
       #	markers_lin.spiral(ax, 'blue', x0, y0, 0.08, 0.08)
       #	
       #for x0,y0 in zip(X[simplemorphtype=='Sp'], log_r_e_eq[simplemorphtype=='Sp']):
       #	markers_lin.spiral(ax, 'blue', x0, y0, 0.08, 0.08)
       #
       #ax.scatter(X[simplemorphtype=='merger'], log_r_e_eq[simplemorphtype=='merger'], marker=r'$\star$', s=300, color='k', **scatter_kwargs)  
	
	ax.scatter(X, log_r_e_eq)
	 
        #plt.xlabel(r'$\alpha \log(\sigma~{\rm [km~s^{-1}]}) + \beta \log(R_{\rm e}~{\rm [kpc]})$', labelpad=13)
        #plt.ylabel(r'$MAG_{\rm sph} \rm ~[mag]$', labelpad=13)
	plt.subplots_adjust(left=0.17,bottom=0.17,right=0.97,top=0.9)
        #plt.axis([-13.5, -5.01, -29, -19])
	
	#ax.text(-13, -20, r'\alpha = ' + )
	
        plt.show()
	#plt.savefig('MAGsph_sigma_Re.pdf', format='pdf', dpi=1000)

	linregr2varEarlyFile.close()
	
def mbh_sigma_resid_Req():
	linregr1varAllFile = open(linregr1varAllFileName)
	linregr1varAll = linregr1varAllFile.readlines()
	i = linregr1varAll.index('log_mbh - log_sigma\n')
	log_sigma_mean = float((linregr1varAll[i+4]).split()[2])
	log_mbh_mean = float((linregr1varAll[i+5]).split()[2])
	A = float((linregr1varAll[i+6]).split()[2])
	B = float((linregr1varAll[i+7]).split()[2])
	
	Y = log_mbh_mean + A + B*(log_sigma-log_sigma_mean)  
	resid = log_mbh - Y
	
        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='E'], resid[simplemorphtype=='E']):
        	markers_lin.elliptical(ax, 'red', (x0), (y0), 0.04, 0.04)
       
        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='E/S0'], resid[simplemorphtype=='E/S0']):
        	markers_lin.lenticular(ax, 'red', x0, y0, 0.07, 0.05)
       
        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='S0'], resid[simplemorphtype=='S0']):
        	markers_lin.lenticular(ax, 'darkorange', x0, y0, 0.07, 0.05)
       
        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='S0/Sp'], resid[simplemorphtype=='S0/Sp']):
        	markers_lin.spiral(ax, 'darkorange', x0, y0, 0.04, 0.04)
        		
        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='Sp'], resid[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.04, 0.04)
        	
        for x0,y0 in zip(log_r_e_eq[simplemorphtype=='Sp'], resid[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.04, 0.04)
       
        ax.scatter(log_r_e_eq[simplemorphtype=='merger'], resid[simplemorphtype=='merger'], marker=r'$\star$', s=300, color='k', **scatter_kwargs)	
 	
	ax.plot([-2,2], [0,0], color='k', ls='--', lw=3)
	
        plt.xlabel(r'$\log(R_{\rm e,eq}~{\rm [kpc]})$', labelpad=13)
        plt.ylabel(r'$M_{\rm BH} - \sigma$ residuals $\rm [dex]$', labelpad=13)
	plt.subplots_adjust(left=0.17,bottom=0.17,right=0.97,top=0.9)
        plt.axis([-1.49, 1.99, -1.49, 1.49])
	
        plt.show()
	#plt.savefig('mbh_sigma_resid_Req.pdf', format='pdf', dpi=1000)

	linregr1varAllFile.close()
	
def mbh_sigma_resid_n_maj():
	linregr1varAllFile = open(linregr1varAllFileName)
	linregr1varAll = linregr1varAllFile.readlines()
	i = linregr1varAll.index('log_mbh - log_sigma\n')
	log_sigma_mean = float((linregr1varAll[i+4]).split()[2])
	log_mbh_mean = float((linregr1varAll[i+5]).split()[2])
	A = float((linregr1varAll[i+6]).split()[2])
	B = float((linregr1varAll[i+7]).split()[2])
	
	Y = log_mbh_mean + A + B*(log_sigma-log_sigma_mean)  
	resid = log_mbh - Y
	
        fig, ax = plt.subplots()

	scatter_kwargs = {"zorder":100}
	error_kwargs = {"lw":.5, "zorder":0}

        for x0,y0 in zip(log_n_maj[simplemorphtype=='E'], resid[simplemorphtype=='E']):
        	markers_lin.elliptical(ax, 'red', (x0), (y0), 0.04, 0.04)
       
        for x0,y0 in zip(log_n_maj[simplemorphtype=='E/S0'], resid[simplemorphtype=='E/S0']):
        	markers_lin.lenticular(ax, 'red', x0, y0, 0.07, 0.05)
       
        for x0,y0 in zip(log_n_maj[simplemorphtype=='S0'], resid[simplemorphtype=='S0']):
        	markers_lin.lenticular(ax, 'darkorange', x0, y0, 0.07, 0.05)
       
        for x0,y0 in zip(log_n_maj[simplemorphtype=='S0/Sp'], resid[simplemorphtype=='S0/Sp']):
        	markers_lin.spiral(ax, 'darkorange', x0, y0, 0.04, 0.04)
        		
        for x0,y0 in zip(log_n_maj[simplemorphtype=='Sp'], resid[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.04, 0.04)
        	
        for x0,y0 in zip(log_n_maj[simplemorphtype=='Sp'], resid[simplemorphtype=='Sp']):
        	markers_lin.spiral(ax, 'blue', x0, y0, 0.04, 0.04)
       
        ax.scatter(log_n_maj[simplemorphtype=='merger'], resid[simplemorphtype=='merger'], marker=r'$\star$', s=300, color='k', **scatter_kwargs)	
 	
	ax.plot([-2,2], [0,0], color='k', ls='--', lw=3)
	
        plt.xlabel(r'$\log(n_{\rm maj}~{\rm [kpc]})$', labelpad=13)
        plt.ylabel(r'$M_{\rm BH} - \sigma$ residuals $\rm [dex]$', labelpad=13)
	plt.subplots_adjust(left=0.17,bottom=0.17,right=0.97,top=0.9)
        plt.axis([-0.49, 1.2, -1.49, 1.49])
	
        plt.show()
	#plt.savefig('mbh_sigma_resid_n_maj.pdf', format='pdf', dpi=1000)

	linregr1varAllFile.close()
	
	

	
#MAGsph_sigma_Req()
#mbh_sigma_resid_Req()
#mbh_sigma_resid_n_maj()
Req_sigma_Ie()

