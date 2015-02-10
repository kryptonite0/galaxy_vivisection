import numpy as np
import sys
from scipy.special import gamma
import os.path

directoriesFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/ProFitErole1.0.1/results/directories.list'

directoriesFile = open(directoriesFileName)
      
directoriesList = []

for line in directoriesFile:
        if (line[0] != '#'):
                directoriesList.append(line.split()[0])

directoriesFile.close()

axisList = ['Rmaj', 'Req']
psfsList = ['gaussian', 'moffat']
scalesList = ['combscale', 'logscale']

fitResultsFile = open('/Users/gsavorgnan/Dropbox/giulia_e_basta/ProFitErole1.0.1/results/1DfitResults.dat', 'w')

sys.stdout = fitResultsFile

print '#', 'galaxy', 
for axis in axisList:
        for psf in psfsList:    
                for scale in scalesList:
                        print 're_sph_' + axis + '_' + psf + '_' + scale,
                        print 'mue_sph_' + axis + '_' + psf + '_' + scale,
                        print 'n_sph_' + axis + '_' + psf + '_' + scale,
                        print 'b_sph_' + axis + '_' + psf + '_' + scale,
                        print 'mtot_sph_' + axis + '_' + psf + '_' + scale,
			print 'h_disc_' + axis + '_' + psf + '_' + scale,
                        print 'mu0_disc_' + axis + '_' + psf + '_' + scale,
                        print 'mtot_disc_' + axis + '_' + psf + '_' + scale,
			print 'Delta_' + axis + '_' + psf + '_' + scale,
print			

for directory in directoriesList:
	print directory,
	for axis in axisList:
		for psf in psfsList:
			for scale in scalesList:
				print '-99  -99  -99  -99  -99  -99  -99  -99  -99 ',
	print
	
        
