import math
from scipy.ndimage.filters import correlate1d
import numpy as np

def gaussFilter(x, y, sigma):

	x = np.asarray(x)
	y = np.asarray(y)
	
	y_conv = [0.0] * (len(x))
	
	for i in range(0,len(x)):
		x0 = x[i]
		kernel = np.exp(-(x - x0)**2/(2*sigma**2))
		norm = sum(kernel)
		kernel_norm = kernel / norm
		
		for j in range(0,len(x)):
			y_conv[i] = y_conv[i] + y[j] * kernel_norm[j]
	
	y_conv = np.asarray(y_conv)
	return y_conv

def moffatFilter(x, y, alpha, beta):

	aa = float(alpha)
	bb = float(beta)

	x = np.asarray(x)
	y = np.asarray(y)

	y_conv = [0.0] * (len(x))
	
	for i in range(0,len(x)):
		x0 = x[i]
		kernel = (1 + ((x - x0)/aa)**2)**(-bb)
		norm = sum(kernel)
		kernel_norm = kernel / norm
		
		for j in range(0,len(x)):
			y_conv[i] = y_conv[i] + y[j] * kernel_norm[j]
	
	y_conv = np.asarray(y_conv)
	return y_conv
	
	
	


def OLDgaussFilter(input, sigma, filterLength, axis = -1, order = 0, output = None,
                      mode = "reflect", cval = 0.0):
    
    """One-dimensional Gaussian filter.

    Parameters
    ----------
    %(input)s
    sigma : scalar
        standard deviation for Gaussian kernel
    filterLength : int
        Length of the filter in units of sigma.	
    %(axis)s
    order : {0, 1, 2, 3}, optional
        An order of 0 corresponds to convolution with a Gaussian
        kernel. An order of 1, 2, or 3 corresponds to convolution with
        the first, second or third derivatives of a Gaussian. Higher
        order derivatives are not implemented
    %(output)s
    %(mode)s
    %(cval)s

    Returns
    -------
    gaussFilter : ndarray

    """
    if order not in range(4):
        raise ValueError('Order outside 0..3 not implemented')
    sd = float(sigma)
    
    
    
    
    
    # make the length of the filter equal to filterLength times the standard
    # deviations:

    lw = int(filterLength * sd + 0.5)  
    
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = 1.0
    sum = 1.0
    sd = sd * sd
    # calculate the kernel:
    for ii in range(1, lw + 1):
        tmp = math.exp(-0.5 * float(ii * ii) / sd)
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum
    # implement first, second and third order derivatives:
    if order == 1 : # first derivative
        weights[lw] = 0.0
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = -x / sd * weights[lw + ii]
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    elif order == 2: # second derivative
        weights[lw] *= -1.0 / sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (x * x / sd - 1.0) * weights[lw + ii] / sd
            weights[lw + ii] = tmp
            weights[lw - ii] = tmp
    elif order == 3: # third derivative
        weights[lw] = 0.0
        sd2 = sd * sd
        for ii in range(1, lw + 1):
            x = float(ii)
            tmp = (3.0 - x * x / sd) * x * weights[lw + ii] / sd2
            weights[lw + ii] = -tmp
            weights[lw - ii] = tmp
    return correlate1d(input, weights, axis, output, mode, cval, 0)

def OLDmoffatFilter(input, alpha, beta, filterLength, axis = -1, output = None,
                      mode = "reflect", cval = 0.0):
    """One-dimensional Moffat filter.

    Parameters
    ----------
    %(input)s
    alpha : scalar
        alpha parameter for Moffat kernel
    beta : scalar
        beta parameter for Moffat kernel	
    %(axis)s
    %(output)s
    %(mode)s
    %(cval)s

    Returns
    -------
    moffatFilter : ndarray

    """
    aa = float(alpha)
    bb = float(beta)
    fwhm = 2*alpha*(2**(1/beta) - 1)**0.5
    # make the length of the filter equal to filterLength times the standard
    # deviations:    
    lw = int(filterLength * fwhm + 0.5)  
    
    weights = [0.0] * (2 * lw + 1)
    weights[lw] = 1.0
    sum = 1.0
    # calculate the kernel:
    for ii in range(1, lw + 1):
        tmp = (1 + (ii/aa)**2)**(-bb)
        weights[lw + ii] = tmp
        weights[lw - ii] = tmp
        sum += 2.0 * tmp
    for ii in range(2 * lw + 1):
        weights[ii] /= sum

    return correlate1d(input, weights, axis, output, mode, cval, 0)

