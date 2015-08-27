from convolution.psfConvolution import doPsfConvolution
from instruments import b_n
import numpy as np

############# BUILD TRUNCATED DISC MODEL ##########

def buildTruncSersic(x, params, component, psfFunction, convolve, zeropoint):  

	namepar1 = 'par1_' + str(component.number) # r_e
	namepar2 = 'par2_' + str(component.number) # mu_e
	namepar3 = 'par3_' + str(component.number) # n
	namepar4 = 'par4_' + str(component.number) # r_out
	
	x_i = x[x<params[namepar4].value]
	x_e = x[x>=params[namepar4].value]
	
	y_tsersic_sb_i = (params[namepar2].value) + (2.5*(b_n.computeb_n(params[namepar3].value))) * ((x_i/(params[namepar1].value))**(1/(params[namepar3].value))-1) / (np.log(10))	
	y_tsersic_i = 10**(0.4*(zeropoint - y_tsersic_sb_i))
	y_tsersic_e = x_e * 0.0
	
	y_tsersic = np.concatenate((y_tsersic_i, y_tsersic_e))
	
	if (convolve):
		y_tsersic = doPsfConvolution(y_tsersic, psfFunction, x)

	return y_tsersic

