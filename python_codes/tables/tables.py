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

path_tables = '/Users/gsavorgnan/galaxy_vivisection/results/tables/'

def make_alister_table1():
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, \
		physres.log_r_e_eq_moffat_comb, errV.perr_log_r_e, errV.merr_log_r_e, \
		physres.mu_e_eq_moffat_comb, errV.perr_mu_e, errV.merr_mu_e, \
		physres.log_n_eq_moffat_comb, errV.perr_log_n, errV.merr_log_n, \
		anc.distance \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	log_r_e = data[1].astype(np.float)
	perr_log_r_e = data[2].astype(np.float)
	merr_log_r_e = data[3].astype(np.float)
	mu_e = data[4].astype(np.float)	
	perr_mu_e = data[5].astype(np.float)		
	merr_mu_e = data[6].astype(np.float)		
	log_n = data[7].astype(np.float)			
	perr_log_n = data[8].astype(np.float)	
	merr_log_n = data[9].astype(np.float)	
	distance = data[10].astype(np.float)	
	
        fName = path_tables + 'r_e-mu_e-n.dat'
        f = open(fName, 'w')
        f.write('# galaxy   distance [Mpc]   log_r_e [log(kpc)]   perr_log_r_e   merr_log_r_e   mu_e [mag/arcsec^2]   perr_mu_e   merr_mu_e   log_n   perr_log_n   merr_log_n \n')
        for g, d, r, pr, mr, m, pm, mm, n, pn, mn in zip(gal_id, distance, log_r_e, perr_log_r_e, merr_log_r_e, mu_e, perr_mu_e, merr_mu_e, log_n, perr_log_n, merr_log_n):
        	f.write(str(g) + '     ' + str(d) + '     ' + str(r) + '     ' + str(pr) + '     ' + str(mr) + '     ' + str(m) + '     ' + str(pm) + '     ' + str(mm) + '     ' + str(n) + '     ' + str(pn) + '     ' + str(mn) + '\n')
        f.close()
       
        ############################
	
	
make_alister_table1()	
	
