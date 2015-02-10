from convolution import filters

######## PERFORM CONVOLUTION WITH GAUSS O moffat ###

def doPsfConvolution(y_component, psfFunction, x):

	if (psfFunction.name == 'gaussian'):
		sigma = psfFunction.gaussianFWHM/2.3548
		#y_component_convolved = filters.gaussFilter(y_component, sigma, filterLength=20, order=0, mode='reflect')
		y_component_convolved = filters.gaussFilter(x, y_component, sigma)
		#y_component_convolved = filters.halfGaussFilter(y_component, sigma, filterLength=40, order=0, mode='reflect')
	elif (psfFunction.name == 'moffat'):
		alpha = psfFunction.moffatAlpha
		beta = psfFunction.moffatBeta
		#y_component_convolved = filters.moffatFilter(y_component, alpha, beta, filterLength=20, mode='reflect')
		y_component_convolved = filters.moffatFilter(x, y_component, alpha, beta)
	return y_component_convolved	
	      
