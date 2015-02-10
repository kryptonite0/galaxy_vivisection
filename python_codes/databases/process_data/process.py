import sqlite3 as sql3
import os
import numpy as np

from conversions import convert 


dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
directoriesFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/1D/galaxies/directories.list'

axisList = ['maj', 'eq']
psfsList = ['gaussian', 'moffat']
scalesList = ['comb', 'log']


def process_single_1D(directory, dist_Mpc):

	collection = []
	collection.append(directory)
	
	for axis in axisList:
		for psf in psfsList:
			for scale in scalesList:
				fitResultsFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/1D/galaxies/' + directory + \
					'/' + directory + '_' + axis + '_' + psf + '_' + scale + '_par_SM.dat'	

				log_r_e = None
				mu_e = None
				log_n = None
				b = None
				mag_sph = None
				mag_tot = None
				
				# get distance
				

				try:
					fitResultsFile = open(fitResultsFileName)
					lines = fitResultsFile.readlines()
					for line in lines:
						if (line.split()[0] == '1') and (line.split()[1] == 'sersic'):
							log_r_e = np.log10(convert.arcsecToKpc(float(line.split()[2]),dist_Mpc))
							mu_e = convert.apparentToAbsoluteMagnitude(float(line.split()[3]),dist_Mpc)
							log_n = np.log10(float(line.split()[4]))
							mag_sph = convert.apparentToAbsoluteMagnitude(float(line.split()[6]),dist_Mpc)
						if (line.split()[0] == 'Total_galaxy_magnitude'):
							mag_tot = convert.apparentToAbsoluteMagnitude(float(line.split()[1]),dist_Mpc)
													
					fitResultsFile.close()
					
				except:
					print 'There was no file named', fitResultsFileName,'.'
	
				collection.append(log_r_e)
				collection.append(mu_e)
				collection.append(log_n)
				collection.append(mag_sph)
				collection.append(mag_tot)
	
	return collection
	
	
def process_all_1D():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	try:
		cur.execute('DROP TABLE OneDFitResultsPhysicalUnits')
		print 'Table OneDFitResultsPhysicalUnits has been erased.'
	except:
		print 'Table OneDFitResultsPhysicalUnits does not exist.'

	directoriesFile = open(directoriesFileName)
	directoriesList = []
	for line in directoriesFile:
		if (line[0] != '#'):
			directoriesList.append(line.split()[0])	
	directoriesFile.close()
	
	create_table_string = 'CREATE TABLE OneDFitResultsPhysicalUnits (gal_id text, ' 
	
	for axis in axisList:
		for psf in psfsList:    
			for scale in scalesList:
				create_table_string += 'log_r_e_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mu_e_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'log_n_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mag_sph_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mag_tot_' + axis + '_' + psf + '_' + scale + ' real, '
	create_table_string = create_table_string[:-2]
	create_table_string += ')'
	
	cur.execute(create_table_string)
	
	for directory in directoriesList:
	
		#get distance in mpc
		cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + directory + '";')
		dist_Mpc = cur.fetchall()[0][0]
		
		collection = process_single_1D(directory, dist_Mpc)
		question_marks = ''
		for i in range(len(collection)-1):
			question_marks += '?,'
		question_marks += '?'	
			
		cur.execute('INSERT INTO OneDFitResultsPhysicalUnits VALUES (' + question_marks + ')', collection)
		
	#cur.execute("select gal_id from OneDFitResults;")
	#results = cur.fetchall()
	#for r in results:
	#	print r
	
	
	connection.commit()
	cur.close()
	connection.close()
	
	print 'Table OneDFitResultsPhysicalUnits has been created.'

	
def main():
	process_all_1D()

main()		
	


