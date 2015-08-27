from convolution.psfConvolution import doPsfConvolution
import numpy as np

############# BUILD TRUNCATED DISC MODEL ##########

def buildTruncDisc(x, params, component, psfFunction, convolve, zeropoint):  

	namepar1 = 'par1_' + str(component.number) # h
	namepar2 = 'par2_' + str(component.number) # mu_0
	namepar3 = 'par3_' + str(component.number) # r_out
	
	x_i = x[x<params[namepar3].value]
	x_e = x[x>=params[namepar3].value]
	
	y_tdisc_sb_i = (params[namepar2].value) + 1.085736*x_i/(params[namepar1].value)
	y_tdisc_i = 10**(0.4*(zeropoint - y_tdisc_sb_i))
	y_tdisc_e = x_e * 0.0
	
	y_tdisc = np.concatenate((y_tdisc_i, y_tdisc_e))
	
	if (convolve):
		y_tdisc = doPsfConvolution(y_tdisc, psfFunction, x)
		
	return y_tdisc

