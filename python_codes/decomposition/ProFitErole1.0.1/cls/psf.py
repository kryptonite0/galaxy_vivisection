from cls import PsfFunction
import numpy as np

def createPsf(psfType,psf_par1,psf_par2,sigmaSmoothing):
    
    psfFunction = None
    
    if (psfType == 'gaussian'): 

        psfFunction = PsfFunction()
        psfFunction.name = 'gaussian'
        psfFunction.gaussianFWHM = psf_par1
        
    if (psfType == 'moffat'): 

        psfFunction = PsfFunction()
        psfFunction.name = 'moffat'
        psfFunction.moffatAlpha = psf_par1  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
        psfFunction.moffatBeta = psf_par2  
        
    gaussianSmoothing = PsfFunction()
    gaussianSmoothing.name = 'gaussian'
    gaussianSmoothing.gaussianFWHM = sigmaSmoothing * 2.3548 
    #gaussianSmoothing.gaussianFWHM = 60 * 2.3548  # m31
    #gaussianSmoothing.gaussianFWHM = 5 * 2.3548     # m81 n5128
    #gaussianSmoothing.gaussianFWHM = 3 * 2.3548    # n4945
    
    return psfFunction, gaussianSmoothing



