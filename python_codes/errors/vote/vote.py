import sqlite3 as sql3
import os
import numpy as np

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
errors1Dfilename = '/Users/gsavorgnan/galaxy_vivisection/data/errors/Errors1D.top'

sigma_mag_sphFileName = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_mag_sph/sigma_mag_sph.txt'
sigma_log_r_eFileName = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_r_e/sigma_log_r_e.txt'
sigma_log_nFileName = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_n/sigma_log_n.txt'
sigma_mu_eFileName = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_mu_e/sigma_mu_e.txt'
sigma_mu_0FileName = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_mu_0/sigma_mu_0.txt'

sigma_mag_sphFile = open(sigma_mag_sphFileName)
lines = sigma_mag_sphFile.readlines()
psigma_mag_sph = np.asarray([float(lines[1].split()[1]), float(lines[1].split()[2]), float(lines[1].split()[3])])
msigma_mag_sph = np.asarray([float(lines[2].split()[1]), float(lines[2].split()[2]), float(lines[2].split()[3])])
sigma_mag_sphFile.close()

sigma_log_r_eFile = open(sigma_log_r_eFileName)
lines = sigma_log_r_eFile.readlines()
psigma_log_r_e = np.asarray([float(lines[1].split()[1]), float(lines[1].split()[2]), float(lines[1].split()[3])])
msigma_log_r_e = np.asarray([float(lines[2].split()[1]), float(lines[2].split()[2]), float(lines[2].split()[3])])
sigma_log_r_eFile.close()

sigma_log_nFile = open(sigma_log_nFileName)
lines = sigma_log_nFile.readlines()
psigma_log_n = np.asarray([float(lines[1].split()[1]), float(lines[1].split()[2]), float(lines[1].split()[3])])
msigma_log_n = np.asarray([float(lines[2].split()[1]), float(lines[2].split()[2]), float(lines[2].split()[3])])
sigma_log_nFile.close()

sigma_mu_eFile = open(sigma_mu_eFileName)
lines = sigma_mu_eFile.readlines()
psigma_mu_e = np.asarray([float(lines[1].split()[1]), float(lines[1].split()[2]), float(lines[1].split()[3])])
msigma_mu_e = np.asarray([float(lines[2].split()[1]), float(lines[2].split()[2]), float(lines[2].split()[3])])
sigma_mu_eFile.close()

sigma_mu_0File = open(sigma_mu_0FileName)
lines = sigma_mu_0File.readlines()
psigma_mu_0 = np.asarray([float(lines[1].split()[1]), float(lines[1].split()[2]), float(lines[1].split()[3])])
msigma_mu_0 = np.asarray([float(lines[2].split()[1]), float(lines[2].split()[2]), float(lines[2].split()[3])])
sigma_mu_0File.close()


def replaceNull(entry):
	if entry == '""':
		entry = None
	return entry	

def errors1D_vote():
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(errors1Dfilename)
	
	try:
		cur.execute('DROP TABLE ErrorsVote')
		print 'Table ErrorsVote has been erased.'
	except:
		print 'Table ErrorsVote does not exist.'
	

	cur.execute('''CREATE TABLE ErrorsVote
		(gal_id text, err_vote integer, 
		perr_log_n real, merr_log_n real, 
		perr_log_r_e real, merr_log_r_e real,  
		perr_mu_e real, merr_mu_e real, 
		perr_mag_sph real, merr_mag_sph real, 
		perr_mu_0 real, merr_mu_0 real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			err_vote = int(line.split()[5])
			
			perr_log_n, merr_log_n = get_error1D_log_n(err_vote)
			perr_log_r_e, merr_log_r_e = get_error1D_log_r_e(err_vote)
			perr_mu_e, merr_mu_e = get_error1D_mu_e(err_vote)
			perr_mag_sph, merr_mag_sph = get_error1D_mag_sph(err_vote)
			perr_mu_0, merr_mu_0 = get_error1D_mu_0(err_vote)
			#perr_mu_0, merr_mu_0 = derive_err_mu_0(gal_id, err_vote, perr_log_n, merr_log_n, perr_mu_e, merr_mu_e)
			
			collection = [gal_id, err_vote, perr_log_n, merr_log_n, perr_log_r_e, merr_log_r_e, perr_mu_e, merr_mu_e, perr_mag_sph, merr_mag_sph, perr_mu_0, merr_mu_0]
			cur.execute('INSERT INTO ErrorsVote VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', collection)

	connection.commit()
	cur.close()
	connection.close()
	data.close()
	
	print 'Table ErrorsVote has been created.'
	
def get_error1D_log_n(vote):

	if (vote > 0) and (vote < 4):	
		index = vote - 1
		perr_log_n = psigma_log_n[index]
		merr_log_n = msigma_log_n[index]
	else:
		perr_log_n = 9999.	
		merr_log_n = 9999.	
	
	return perr_log_n, merr_log_n

def get_error1D_log_r_e(vote):
	
	if (vote > 0) and (vote < 4):	
		index = vote - 1
		perr_log_r_e = psigma_log_r_e[index]
		merr_log_r_e = msigma_log_r_e[index]
	else:
		perr_log_r_e = 9999.
		merr_log_r_e = 9999.
			
	return perr_log_r_e, merr_log_r_e

def get_error1D_mu_e(vote):
	
	if (vote > 0) and (vote < 4):	
		index = vote - 1
		perr_mu_e = psigma_mu_e[index]
		merr_mu_e = msigma_mu_e[index]
	else:
		perr_mu_e = 9999.	
		merr_mu_e = 9999.	
	
	return perr_mu_e, merr_mu_e

def get_error1D_mag_sph(vote):
	
	if (vote > 0) and (vote < 4):	
		index = vote - 1
		perr_mag_sph = psigma_mag_sph[index]
		merr_mag_sph = msigma_mag_sph[index]
	else:
		perr_mag_sph = 9999.	
		merr_mag_sph = 9999.	
	
	return perr_mag_sph, merr_mag_sph
	
def derive_err_mu_0(gal_id, vote, perr_log_n, merr_log_n, perr_mu_e, merr_mu_e):

	if (vote > 0) and (vote < 4):
		connection = sql3.connect(dbname)
		cur = connection.cursor()
	
		query = 'SELECT one.n_maj_moffat_comb \
			FROM OneDFitResults AS one \
			WHERE one.gal_id = "' + gal_id + '";'
			
		cur.execute(query)
		n = cur.fetchall()[0][0]
	
		perr_n = n*(10**perr_log_n - 1)
		merr_n = n*(1 - 10**(-merr_log_n))
	
		perr_mu_0 = (perr_mu_e**2 + ((2.5*1.9992/np.log(10))**2)*merr_n**2)**0.5
		merr_mu_0 = (merr_mu_e**2 + ((2.5*1.9992/np.log(10))**2)*perr_n**2)**0.5
		#print gal_id, vote, n, perr_n, merr_n, perr_mu_e, merr_mu_e, perr_mu_0,merr_mu_0	

	else:
		perr_mu_0,merr_mu_0 = 9999., 9999.
	
	return perr_mu_0,merr_mu_0
	
def get_error1D_mu_0(vote):
	
	if (vote > 0) and (vote < 4):	
		index = vote - 1
		perr_mu_0 = psigma_mu_0[index]
		merr_mu_0 = msigma_mu_0[index]
	else:
		perr_mu_0 = 9999.	
		merr_mu_0 = 9999.	
	
	return perr_mu_0, merr_mu_0
	
def main():
	errors1D_vote()	
	
main()
