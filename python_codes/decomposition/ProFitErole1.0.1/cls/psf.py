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
        
        psfList = [gaussianPsf, moffatPsf]
#         psfList = [moffatPsf]

#     if (Settings.observation == 'n1277'):
#     
#         moffatPsf = PsfFunction()
#         moffatPsf.name = 'moffat'
#         moffatPsf.moffatAlpha = 2.6/(2*np.sqrt(2**(1/2.9)-1))  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
#         moffatPsf.moffatBeta = 2.9  
#         
#         psfList = [moffatPsf]

    gaussianSmoothing = PsfFunction()
    gaussianSmoothing.name = 'gaussian'
    gaussianSmoothing.gaussianFWHM = 60 * 2.3548 * Settings.pxlToArcsec 
    
    return psfList, gaussianSmoothing



