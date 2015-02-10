from convolution.psfConvolution import doPsfConvolution
import numpy as np


def buildGaussian(x, params, component, psfFunction, convolve, zeropoint):

	namepar1 = 'par1_' + str(component.number) # fwhm
	namepar2 = 'par2_' + str(component.number) # mag
	y_gauss = 10**(zeropoint/2.5) * (10**(-params[namepar2].value/2.5) * np.exp(-x**2/(2*(params[namepar1].value/2.3548)**2)))
	#print y_gauss
	if (convolve):
		y_gauss = doPsfConvolution(y_gauss, psfFunction, x)
	return y_gauss

def computeGaussianParameters(fwhm, mu_0):

	sigma = fwhm/2.3548
	m_tot = mu_0 - 2.5*np.log10(2*np.pi*(sigma)**2)
	
	return m_tot