import numpy as np
import sys
from scipy.special import gamma
import os.path

directoriesFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/directories.list'

directoriesFile = open(directoriesFileName)

directoriesList = []

for line in directoriesFile:
	if (line[0] != '#'):
		directoriesList.append(line.split()[0])

directoriesFile.close()

axisList = ['maj', 'eq']
psfsList = ['gaussian', 'moffat']
#scalesList = ['logscale', 'linscale1px', 'linscale05px']
scalesList = ['comb', 'log']

fitResultsFile = open('/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/1DfitsReloadResults.dat', 'w')

sys.stdout = fitResultsFile

print '#', 'galaxy', 
for axis in axisList:
	for psf in psfsList:    
		for scale in scalesList:
			print 're_' + axis + '_' + psf + '_' + scale,
			print 'mue_' + axis + '_' + psf + '_' + scale,
			print 'n_' + axis + '_' + psf + '_' + scale,
			print 'msph_' + axis + '_' + psf + '_' + scale,
			print 'mtot_' + axis + '_' + psf + '_' + scale,          
print 


for directory in directoriesList:

	print directory,
	
	for axis in axisList:
							
		for psf in psfsList:
			
			for scale in scalesList:
					
				outputFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/' + directory + '/' + directory + '_' + axis + '_' + psf + '_' + scale + '_par_SM.dat'
				
				r_e = '""'
				mu_e = '""'
				n = '""'
				b = '""'
				mag_sph = '""'
				mag_tot = '""'

				try:
					outputFile = open(outputFileName)
					lines = outputFile.readlines()
					for line in lines:
						if (line.split()[0] == '1') and (line.split()[1] == 'sersic'):
							r_e = float(line.split()[2])
							mu_e = float(line.split()[3])
							n = float(line.split()[4])
							mag_sph = float(line.split()[6])
						if (line.split()[0] == 'Total_galaxy_magnitude'):
							mag_tot = float(line.split()[1])
													
					outputfile.close()
					
				except:
					None
				
# 				if not os.path.exists(outputFileName):
				
				print r_e, mu_e, n , mag_sph, mag_tot,
					
	print '\n',					
