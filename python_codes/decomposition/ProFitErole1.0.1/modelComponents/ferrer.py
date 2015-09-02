from convolution.psfConvolution import doPsfConvolution
import numpy as np
from scipy import special #.hyp2f1

def buildFerrer(x, params, component, psfFunction, convolve, zeropoint):

	namepar1 = 'par1_' + str(component.number) # r_out
	namepar2 = 'par2_' + str(component.number) # mu_0
	namepar3 = 'par3_' + str(component.number) # alpha
	namepar4 = 'par4_' + str(component.number) # beta
	
	x_i = x[x<params[namepar1].value]
	x_e = x[x>=params[namepar1].value]
	
	y_ferrer_i = 10**(zeropoint/2.5) * (10**(-params[namepar2].value/2.5)) * ( 1 - (x_i/params[namepar1].value)**(2 - (params[namepar4].value)) )**(params[namepar3].value)
	y_ferrer_e = x_e * 0.0
	y_ferrer = np.concatenate((y_ferrer_i, y_ferrer_e))
			
	if (convolve):
		y_ferrer = doPsfConvolution(y_ferrer, psfFunction, x)
	return y_ferrer

def computeFerrerParameters(r_out, mu_0, alpha, beta):
	
	m_tot = mu_0 - 2.5 * np.log10(2 * np.pi * r_out**2 * (0.5 * special.hyp2f1(-alpha, 2/(2-beta), (-beta)/(2-beta), 1)))
	
	return m_tot
