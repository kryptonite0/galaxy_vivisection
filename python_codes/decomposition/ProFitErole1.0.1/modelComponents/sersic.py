import numpy as np
from instruments import b_n
from convolution.psfConvolution import doPsfConvolution
from scipy.special import gammainc, gamma

def buildSersic(x, params, component, psfFunction, convolve, zeropoint):

	namepar1 = 'par1_' + str(component.number) # r_e
	namepar2 = 'par2_' + str(component.number) # mu_e
	namepar3 = 'par3_' + str(component.number) # n
	
	y_sersic_sb = (params[namepar2].value) + (2.5*(b_n.computeb_n(params[namepar3].value))) * ((x/(params[namepar1].value))**(1/(params[namepar3].value))-1) / (np.log(10))	
	y_sersic = 10**(0.4*(zeropoint - y_sersic_sb))
	
	if (convolve):
		y_sersic = doPsfConvolution(y_sersic, psfFunction, x)

	return y_sersic

def computeSersicParameters(mu_e, r_e, n):
	
	b = b_n.computeb_n(n)
	m_tot = mu_e - 5*np.log10(r_e) - 2.5*np.log10(2*np.pi*n*np.exp(b)*gamma(2*n)/(b**(2*n)))
	
	return b, m_tot
	
def compute_mu_e(m_tot, r_e, n):

	b = b_n.computeb_n(n)
	mu_e = m_tot + 5*np.log10(r_e) + 2.5*np.log10(2*np.pi*n*np.exp(b)*gamma(2*n)/(b**(2*n)))
	
	return mu_e

