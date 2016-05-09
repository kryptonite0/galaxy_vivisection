from cls.cls import PsfFunction
from cls.psf import createPsf
import numpy as np
    
def readInitialSettings(txt, Settings):

    data = open(txt)
    
    lines = data.readlines()
    
    galaxyName = str((lines[0]).split()[0])
    prefixEllipseOutput = str((lines[1]).split()[0])
    skyRMS = float((lines[2]).split()[0])
    Settings.pxlToArcsec = float((lines[3]).split()[0])
    Settings.zeropoint = float((lines[4]).split()[0])
    psfType = str((lines[5]).split()[0])
    psf_par1 = Settings.pxlToArcsec*float((lines[6]).split()[0])
    psf_par2 = float((lines[7]).split()[1])
    smoothing = str((lines[8]).split()[0])
    sigmaSmoothing = float((lines[9]).split()[0])
        
    if smoothing == 'yes':
    	Settings.smoothing == True
	
    psfFunction, gaussianSmoothing = createPsf(Settings,psfType,psf_par1,psf_par2,sigmaSmoothing)

    return galaxyName, prefixEllipseOutput, skyRMS, Settings, psfFunction, gaussianSmoothing	
	
	
