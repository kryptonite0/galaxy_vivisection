from convolution.psfConvolution import doPsfConvolution
import numpy as np
from scipy import special

def buildGaussRing(x, params, component, psfFunction, convolve, zeropoint):
	
	namepar1 = 'par1_' + str(component.number) # fwhm
	namepar2 = 'par2_' + str(component.number) # mag
	namepar3 = 'par3_' + str(component.number) # r_0
	
	y_gring = 10**(zeropoint/2.5) * (10**(-params[namepar2].value/2.5) * np.exp(-(x-params[namepar3].value)**2/(2*(params[namepar1].value/2.3548)**2)))
	#print y_gauss
	if (convolve):
		y_gring = doPsfConvolution(y_gring, psfFunction, x)
	return y_gring

def computeGaussRingParameters(fwhm, mu_0, r_0):
	
	sigma = fwhm/2.3548	
	m_tot = mu_0 - 2.5*np.log10(2**0.5 * np.pi * sigma *(2**0.5 * sigma * np.exp(-r_0**2 / (2*sigma**2)) + r_0 * np.pi**0.5 * (1 - special.erf(-r_0 / (2**0.5 * sigma)))))
	
	return m_tot