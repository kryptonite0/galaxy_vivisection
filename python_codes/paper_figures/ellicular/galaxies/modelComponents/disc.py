from convolution.psfConvolution import doPsfConvolution
import numpy as np
from instruments import b_n
from scipy.special import gammainc, gamma


def buildDisc(x, params, component, psfFunction, convolve, zeropoint):  

	namepar1 = 'par1_' + str(component.number) # h
	namepar2 = 'par2_' + str(component.number) # mu_0
	y_disc_sb = (params[namepar2].value) + 1.085736*x/(params[namepar1].value)
	y_disc = 10**(0.4*(zeropoint - y_disc_sb))
	
	if (convolve):
		y_disc = doPsfConvolution(y_disc, psfFunction, x)
		
	return y_disc

def computeDiscParameters(mu_0, h):
	r_e = 1.678*h
	mu_e = mu_0 + 1.821865
	n = 1
	b=1.678
	m_tot = mu_e - 5*np.log10(r_e) - 2.5*np.log10(2*np.pi*n*np.exp(b)*gamma(2*n)/(b**(2*n)))
	
	return m_tot