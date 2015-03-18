import sqlite3 as sql3
import os
import numpy as np

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
directoriesFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/1D/galaxies/directories.list' # for 1D fits
bestfitmodelFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/2D/galaxies/bestfitmodel.list' # for 2D fits

axisList = ['maj', 'eq']
psfsList = ['gaussian', 'moffat']
scalesList = ['comb', 'log']

zeropoint = 17.25585
pxlToArcsec = 1.2233

def collect_single_1D(directory):

	collection = []
	collection.append(directory)
	
	for axis in axisList:
		for psf in psfsList:
			for scale in scalesList:
				fitResultsFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/1D/galaxies/' + directory + \
					'/' + directory + '_' + axis + '_' + psf + '_' + scale + '_par_SM.dat'	

				r_e = None
				mu_e = None
				n = None
				b = None
				mag_sph = None
				mag_tot = None

				try:
					fitResultsFile = open(fitResultsFileName)
					lines = fitResultsFile.readlines()
					for line in lines:
						if (line.split()[0] == '1') and (line.split()[1] == 'sersic'):
							r_e = float(line.split()[2])
							mu_e = float(line.split()[3])
							n = float(line.split()[4])
							mag_sph = float(line.split()[6])
						if (line.split()[0] == 'Total_galaxy_magnitude'):
							mag_tot = float(line.split()[1])
						if (line.split()[0] == 'Delta'):
							delta = float(line.split()[1])	
													
					fitResultsFile.close()
					
				except:
					print 'There was no file named', fitResultsFileName,'.'
	
				collection.append(r_e)
				collection.append(mu_e)
				collection.append(n)
				collection.append(mag_sph)
				collection.append(mag_tot)
				collection.append(delta)
	
	return collection
	
def collect_single_2D(directory, bestfitFile, mosaic):

	collection = []
	collection.append(directory)
	
	for psf in psfsList:
		fitParametersFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/2D/galaxies/' + directory + '/' + \
			bestfitFile + '_' + psf + '_' + mosaic
		fitFluxesFileName = '/Users/gsavorgnan/galaxy_vivisection/analysis/2D/galaxies/' + directory + '/totalFluxes_' + \
			bestfitFile + '_' + psf + '_' + mosaic	
		
		if os.path.exists(fitParametersFileName):
			fitParametersFile = open(fitParametersFileName)
			lines = fitParametersFile.readlines()
			
			#X0 = float(lines[4].split()[1])
			#Y0 = float(lines[5].split()[1])
			#PA = float(lines[7].split()[1])
			#ell = float(lines[8].split()[1])
			#c0 = float(lines[9].split()[1]) 
			n = float(lines[10].split()[1])
			I_e = float(lines[11].split()[1])
			r_e_pix = float(lines[12].split()[1])
			
			mu_e = zeropoint - 2.5 * np.log10(I_e)
			r_e = r_e_pix*pxlToArcsec
			
			fitParametersFile.close()
			
		else:
			#X0 = 
			#Y0 = 
			#PA = 
			#ell = 
			#c0 = 
			n = None 
			mu_e = None 
			r_e = None 
			print 'There was no file named', fitParametersFileName,'.'

		collection.append(r_e)
		collection.append(mu_e)
		collection.append(n)
		
		mag_sph = None
		if os.path.exists(fitFluxesFileName):
			fitFluxesFile = open(fitFluxesFileName)
			lines = fitFluxesFile.readlines()
			
			for line in lines:
				
				if 'Sersic' in line and 'Function' not in line:
					#print line
					flux_sph = float(line.split()[1])
					mag_sph = float(line.split()[2])
					break
				#if ('Exponential' in line or 'EdgeOnDisk' in line) and 'Function' not in line:
				#	#print line
				#	#flux_disc = line.split()[1]
				#	mag_disc = line.split()[2]
				#if 'FlatSky' in line and 'Function' not in line:
				#	#print line
				#	flux_sky = float(line.split()[1])
				#	#mag_sky = line.split()[2]
				#if 'Total' in line:
				#	flux_total = float(line.split()[1])
			
			#if flux_total>0 and flux_sky>0:
			#	flux_galaxy = flux_total - flux_sky
			#	fraction_sph = flux_sph/flux_galaxy
			fitFluxesFile.close()	
					
		else:
			print 'There was no file named', fitFluxesFileName,'.'
				
		collection.append(mag_sph)
	
	return collection
	
	
def collect_all_1D():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	try:
		cur.execute('DROP TABLE OneDFitResults')
		print 'Table OneDFitResults has been erased.'
	except:
		print 'Table OneDFitResults does not exist.'

	#directoriesFile = open(directoriesFileName)
	#directoriesList = []
	#for line in directoriesFile:
	#	if (line[0] != '#'):
	#		directoriesList.append(line.split()[0])	
	#directoriesFile.close()
	
	create_table_string = 'CREATE TABLE OneDFitResults (gal_id text, ' 
	
	for axis in axisList:
		for psf in psfsList:    
			for scale in scalesList:
				create_table_string += 'r_e_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mu_e_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'n_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mag_sph_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'mag_tot_' + axis + '_' + psf + '_' + scale + ' real, '
				create_table_string += 'delta_' + axis + '_' + psf + '_' + scale + ' real, '
	create_table_string = create_table_string[:-2]
	create_table_string += ')'
	
	cur.execute(create_table_string)
	
	cur.execute('SELECT gal_id FROM Ancillary WHERE fit1D_done = 1;')
	directoriesList = cur.fetchall()
	#print 'dir', directoriesList
	
	for directory in directoriesList:
		
		collection = collect_single_1D(directory[0])
		question_marks = ''
		for i in range(len(collection)-1):
			question_marks += '?,'
		question_marks += '?'	
			
		cur.execute('INSERT INTO OneDFitResults VALUES (' + question_marks + ')', collection)
		
	#cur.execute("select gal_id from OneDFitResults;")
	#results = cur.fetchall()
	#for r in results:
	#	print r
	
	
	connection.commit()
	cur.close()
	connection.close()
	
	print 'Table OneDFitResults has been created.'

def collect_all_2D():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	try:
		cur.execute('DROP TABLE TwoDFitResults')
		print 'Table TwoDFitResults has been erased.'
	except:
		print 'Table TwoDFitResults does not exist.'

	bestfitmodelFile = open(bestfitmodelFileName)
	directoriesList = []
	bestfitFilesList = []
	mosaicsList = []
	
	for line in bestfitmodelFile:
	        if (line[0] != '#'):
	                directoriesList.append(line.split()[0])
			bestfitFilesList.append(line.split()[1])
			mosaicsList.append(line.split()[2])
	bestfitmodelFile.close()
	
	create_table_string = 'CREATE TABLE TwoDFitResults (gal_id text, r_e_gaussian real, mu_e_gaussian real, \
		n_gaussian real, mag_sph_gaussian real, \
		r_e_moffat real, mu_e_moffat real, n_moffat real, mag_sph_moffat real)' 
	
	cur.execute(create_table_string)
	
	for i in range(0,len(directoriesList)):
		directory = directoriesList[i]
		bestfitFile = bestfitFilesList[i]
		mosaic = mosaicsList[i]
		
		collection = collect_single_2D(directory, bestfitFile, mosaic)
		#print directory, collection
		cur.execute('INSERT INTO TwoDFitResults VALUES (?,?,?,?,?,?,?,?,?)', collection)
		
	#cur.execute("select gal_id from OneDFitResults;")
	#results = cur.fetchall()
	#for r in results:
	#	print r
	
	
	connection.commit()
	cur.close()
	connection.close()
	
	print 'Table TwoDFitResults has been created.'

	
def main():
	collect_all_1D()
	collect_all_2D()

main()		
	
	
	
	
		

