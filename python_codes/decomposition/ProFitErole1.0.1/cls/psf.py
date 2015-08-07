from cls import PsfFunction
import numpy as np

def createPsf(Settings):
    
    if (Settings.observation == 'spitzer3.6um'): 

        gaussianPsf = PsfFunction()
        gaussianPsf.name = 'gaussian'
        gaussianPsf.gaussianFWHM = 1.6648 * Settings.pxlToArcsec 
        
        moffatPsf = PsfFunction()
        moffatPsf.name = 'moffat'
        moffatPsf.moffatAlpha = (1.61467/(2*np.sqrt(2**(1/4.39)-1)) ) * Settings.pxlToArcsec  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
        moffatPsf.moffatBeta = 4.39  
        
#        psfList = [gaussianPsf, moffatPsf]
        psfList = [moffatPsf]

    if (Settings.observation == 'n1271'):
    
        gaussianPsf = PsfFunction()
        gaussianPsf.name = 'gaussian'
        gaussianPsf.gaussianFWHM = 2.3 * Settings.pxlToArcsec 
        
        psfList = [gaussianPsf]

    if (Settings.observation == 'n4342'):
    
        moffatPsf = PsfFunction()
        moffatPsf.name = 'moffat'
        moffatPsf.moffatAlpha = (2./(2*np.sqrt(2**(1/1.9)-1)) ) * Settings.pxlToArcsec  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
        moffatPsf.moffatBeta = 1.9  
        
        psfList = [moffatPsf]

    if (Settings.observation == 'LEDA'):
    
        moffatPsf = PsfFunction()
        moffatPsf.name = 'moffat'
        moffatPsf.moffatAlpha = (2.69/(2*np.sqrt(2**(1/6.62)-1)) ) * Settings.pxlToArcsec  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
        moffatPsf.moffatBeta = 6.62  
        
        psfList = [moffatPsf]

    gaussianSmoothing = PsfFunction()
    gaussianSmoothing.name = 'gaussian'
    #gaussianSmoothing.gaussianFWHM = 60 * 2.3548 * Settings.pxlToArcsec  # m31
    gaussianSmoothing.gaussianFWHM = 5 * 2.3548 * Settings.pxlToArcsec    # m81 n5128
    #gaussianSmoothing.gaussianFWHM = 3 * 2.3548 * Settings.pxlToArcsec    # n4945
    
    return psfList, gaussianSmoothing



