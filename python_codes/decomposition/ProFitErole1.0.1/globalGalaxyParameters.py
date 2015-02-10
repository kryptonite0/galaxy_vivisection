import numpy as np
import sys
from scipy.special import gamma
import os.path
from lmfit import Parameters

from reading.bestFitModel import readBestFitModel
from instruments.model import buildModel
from cls.cls import PsfFunction

class Settings:
    
    smoothing = False
    
    pxlToArcsec = 1.2232836
    zeropoint = 17.25585 #spitzer 3.6 um

    minx = 0.
    maxx = 10000.
    stepx = 0.1
    
    moffatPsf = PsfFunction()
    moffatPsf.name = 'moffat'
    moffatPsf.moffatAlpha = (1.61467/(2*np.sqrt(2**(1/4.39)-1)) ) * pxlToArcsec  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
    moffatPsf.moffatBeta = 4.39  


############# main body here ##############

directoriesFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/directories.list'

directoriesFile = open(directoriesFileName)

directoriesList = []

for line in directoriesFile:
    if (line[0] != '#'):
        directoriesList.append(line.split()[0])

directoriesFile.close()

outputFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/globalGalaxyParameters.dat'
outputFile = open(outputFileName, 'w')
sys.stdout = outputFile
print '# galaxy  R_e_tot  mag_tot'

for directory in directoriesList:
    
    bestFitModelFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/1DfitsReload/' + directory + '/' + directory + '_eq_moffat_comb_par_SM.dat'
    
    componentslist, params = readBestFitModel(bestFitModelFileName)
    
#     xxx = np.arange(Settings.minx+Settings.stepx/2,Settings.maxx+Settings.stepx/2,Settings.stepx)
    xxx = np.arange(Settings.minx,Settings.maxx,Settings.stepx)
    emptyList = []
    
    model_sb, good_model_sb = buildModel(params, componentslist, xxx, emptyList, Settings.moffatPsf, False, None, None, Settings, 'comb', None)
    
    model = 10**(-0.4*(model_sb - Settings.zeropoint))
    
#     totalLuminosity = np.sum(model*Settings.stepx*2*np.pi*xxx)
    
    totalLuminosity = 0
    for i in np.arange(0,len(xxx)-1,1):
        totalLuminosity = totalLuminosity + (model[i]+model[i+1])*0.5*Settings.stepx*2*np.pi*xxx[i]
    
    halfLuminosity = totalLuminosity*0.5

    tmp = 0
    for i in np.arange(0,len(xxx)-1,1):
        if (tmp < halfLuminosity):
            tmp = tmp + (model[i]+model[i+1])*0.5*Settings.stepx*2*np.pi*xxx[i]
        else:
            R_e = xxx[i]
            break
    
    totalMagnitude = Settings.zeropoint - 2.5*np.log10(totalLuminosity)
    
    
    print directory, R_e, totalMagnitude            

# 
# build tot model on (0,+inf)
#     build single functions 
#     sum them together
# 
# integrate tot light
# 
# integrate up to 0.5 tot light to get R_e of galaxy