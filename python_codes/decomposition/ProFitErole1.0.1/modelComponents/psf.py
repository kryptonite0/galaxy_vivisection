import numpy as np

def buildPsf(x, params, component, psfFunction, zeropoint):

	if (psfFunction.name == 'gaussian'):
		namepar2 = 'par2_' + str(component.number) # mu_0
		fwhm = psfFunction.gaussianFWHM
		y_psf = 10**(zeropoint/2.5) * (10**(-params[namepar2].value/2.5) * np.exp(-x**2/(2*(fwhm/2.3548)**2)))
	elif (psfFunction.name == 'moffat'):
		namepar2 = 'par2_' + str(component.number) # mu_0
		alpha = psfFunction.moffatAlpha
		beta = psfFunction.moffatBeta
		y_psf = 10**(zeropoint/2.5) * (10**(-params[namepar2].value/2.5) * (1 + (x/alpha)**2)**(-beta))
	
	return y_psf	

