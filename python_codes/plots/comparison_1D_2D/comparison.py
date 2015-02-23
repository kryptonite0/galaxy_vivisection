import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 32})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
path = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_1D_2D/'
paper_path = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/images/'

def comparison_r_e():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, \
		one.r_e_maj_moffat_comb, \
		two.r_e_moffat \
		FROM Ancillary as anc \
		JOIN OneDFitResults as one ON anc.gal_id = one.gal_id\
		JOIN TwoDFitResults AS two ON anc.gal_id = two.gal_id \
		WHERE anc.fit1D_done = 1 AND anc.fit2D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	r_e_maj_1D = data[1].astype(np.float)
	r_e_2D = data[2].astype(np.float)
	
	#print r_e_maj_1D, r_e_2D
	print 'there are ', len(r_e_maj_1D), ' points in this plot'
	
	fig, ax = plt.subplots()
	ax.scatter(r_e_maj_1D, r_e_2D, s=120, c='black')
	ax.set_xscale('log')
	ax.set_yscale('log')
        ticks = (np.asarray([1., 10., 100.]))
        ticks_labels = ['$1$','$10$','$100$']
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks_labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks_labels)
	ax.axis([0.5,700.,0.5,700.])
	ax.plot([0.1,1000], [0.1,1000], 'k--', linewidth=3)
        plt.xlabel(r'$R_{\rm e,sph}^{\rm 1D,maj} \rm~[arcsec]$', labelpad=10)
        plt.ylabel(r'$R_{\rm e,sph}^{\rm 2D,maj} \rm~[arcsec]$', labelpad=10)
	plt.subplots_adjust(left=0.2,bottom=0.2,right=0.98,top=0.98)

	#plt.show()
	plt.savefig(path + 'comparison_r_e.eps', format='eps', dpi=1000)
	plt.savefig(paper_path + 'comparison_r_e.eps', format='eps', dpi=1000)
	
def comparison_n():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, \
		one.n_eq_moffat_comb, \
		two.n_moffat \
		FROM Ancillary as anc \
		JOIN OneDFitResults as one ON anc.gal_id = one.gal_id\
		JOIN TwoDFitResults AS two ON anc.gal_id = two.gal_id \
		WHERE anc.fit1D_done = 1 AND anc.fit2D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	n_eq_1D = data[1].astype(np.float)
	n_2D = data[2].astype(np.float)
	
	fig, ax = plt.subplots()
	ax.scatter(n_eq_1D, n_2D, s=120, c='black')
	ax.set_xscale('log')
	ax.set_yscale('log')
        ticks = (np.asarray([1., 10.]))
        ticks_labels = ['$1$','$10$']
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks_labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks_labels)
	ax.axis([0.7,12.,0.7,12.])
	ax.plot([0.1,100], [0.1,100], 'k--', linewidth=3)
        plt.xlabel(r'$n_{\rm sph}^{\rm 1D,eq}$', labelpad=10)
        plt.ylabel(r'$n_{\rm sph}^{\rm 2D}$', labelpad=10)
	plt.subplots_adjust(left=0.15,bottom=0.2,right=0.98,top=0.98)

	plt.savefig(path + 'comparison_n.eps', format='eps', dpi=1000)
	plt.savefig(paper_path + 'comparison_n.eps', format='eps', dpi=1000)
	#plt.show()

def comparison_mag_sph():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, \
		one.mag_sph_eq_moffat_comb, \
		two.mag_sph_moffat \
		FROM Ancillary as anc \
		JOIN OneDFitResults as one ON anc.gal_id = one.gal_id\
		JOIN TwoDFitResults AS two ON anc.gal_id = two.gal_id \
		WHERE anc.fit1D_done = 1 AND anc.fit2D_done = 1;'
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	mag_sph_eq_1D = data[1].astype(np.float)
	mag_sph_2D = data[2].astype(np.float)
	
	fig, ax = plt.subplots()
	ax.scatter(mag_sph_eq_1D, mag_sph_2D, s=120, c='black')
	#ax.set_xscale('log')
	#ax.set_yscale('log')
        #ticks = (np.asarray([1., 10.]))
        #ticks_labels = ['$1$','$10$']
        #ax.set_xticks(ticks)
        #ax.set_xticklabels(ticks_labels)
        #ax.set_yticks(ticks)
        #ax.set_yticklabels(ticks_labels)
	ax.axis([4.1,10.9,4.1,10.9])
	ax.plot([0.1,100], [0.1,100], 'k--', linewidth=3)
        plt.xlabel(r'$m_{\rm sph}^{\rm 1D} \rm~[mag]$', labelpad=10)
        plt.ylabel(r'$m_{\rm sph}^{\rm 2D} \rm~[mag]$', labelpad=10)
	plt.subplots_adjust(left=0.15,bottom=0.2,right=0.98,top=0.98)

	plt.savefig(path + 'comparison_mag.eps', format='eps', dpi=1000)
	plt.savefig(paper_path + 'comparison_mag.eps', format='eps', dpi=1000)
	#plt.show()


def main():
	comparison_r_e()
	comparison_n()
	comparison_mag_sph()
	
main()		
