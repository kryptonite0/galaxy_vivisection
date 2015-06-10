import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from scipy import stats
#from instruments.linear_regression import bces
#from instruments.linear_regression import predband
#from instruments.linear_regression import fitexy
from instruments import b_n
#from instruments.linear_regression import colorline


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 12})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

# get data
connection = sql3.connect(dbname)
cur = connection.cursor()

getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.bar, \
	anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.sigma, \
	physres.mag_sph_eq_moffat_comb, physres.mag_tot_eq_moffat_comb, \
	physres.log_n_maj_moffat_comb, physres.log_n_eq_moffat_comb, \
	physres.mu_e_maj_moffat_comb, physres.mu_e_eq_moffat_comb, \
	physres.log_r_e_maj_moffat_comb, physres.log_r_e_eq_moffat_comb \
	FROM Ancillary AS anc \
	JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
	WHERE anc.fit1D_done = 1;'

cur.execute(getdata_query)
datalist = cur.fetchall()
data= np.asarray(datalist).transpose()
#print data
gal_id = data[0]
simplemorphtype = data[1]
core = data[2].astype(np.int)
bar = data[3].astype(np.int)

mbh = data[4].astype(np.float)
log_mbh = np.log10(mbh)

sigma = data[7].astype(np.float)
# assign value to n3079
sigma[gal_id=='n3079'] = 105
log_sigma = np.log10(sigma)

mag_sph = data[8].astype(np.float)
mag_tot = data[9].astype(np.float)

log_n_maj = data[10].astype(np.float)
n_maj = 10**log_n_maj
log_n_eq = data[11].astype(np.float)
n_eq = 10**log_n_eq

mu_e_maj = data[12].astype(np.float)
mu_e_eq = data[13].astype(np.float)

# compute mu_0
b_maj = mu_e_maj * [0.0]
for i in range(len(b_maj)):
	b_maj[i] = b_n.computeb_n(n_maj[i])
b_eq = mu_e_eq * [0.0]
for i in range(len(b_eq)):
	b_eq[i] = b_n.computeb_n(n_eq[i])

mu_0_maj = mu_e_maj - 2.5*b_maj/np.log(10)
mu_0_eq = mu_e_eq - 2.5*b_eq/np.log(10)

log_r_e_maj = data[14].astype(np.float)
r_e_maj = 10**log_r_e_maj
log_r_e_eq = data[15].astype(np.float)
r_e_eq = 10**log_r_e_eq

# compute mean vector
mean_vector = np.array([[np.mean(log_mbh)], [np.mean(log_sigma)], [np.mean(mag_sph)], [np.mean(mag_tot)], 
	[np.mean(log_n_maj)], [np.mean(log_n_eq)], [np.mean(mu_e_maj)], [np.mean(mu_e_eq)], [np.mean(mu_0_maj)], 
	[np.mean(mu_0_eq)], [np.mean(log_r_e_maj)], [np.mean(log_r_e_eq)]])
print('Mean Vector:\n', mean_vector)

# compute covariance matrix
cov_mat = np.cov([log_mbh, log_sigma, mag_sph, mag_tot, log_n_maj, log_n_eq, mu_e_maj, mu_e_eq, mu_0_maj, mu_0_eq, log_r_e_maj, log_r_e_eq])
print('Covariance Matrix:\n', cov_mat)

# eigenvectors and eigenvalues for the from the covariance matrix
eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)
print eig_val_cov, eig_vec_cov


