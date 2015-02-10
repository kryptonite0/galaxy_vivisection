import numpy
import scipy

def scatterfit(x,y,a,b):  
    """  
    Compute the mean deviation of the data about the linear model given if A,B  
    (y=ax+b) provided as arguments. Otherwise, compute the mean deviation about   
    the best-fit line.  
   
    x,y assumed to be Numpy arrays. a,b scalars.  
    Returns the float sd with the mean deviation.  
   
    Author: Rodrigo Nemmen  
    """  
      
    # Std. deviation of an individual measurement (Bevington, eq. 6.15)  
    N=numpy.size(x)  
    sd=1./(N-2.)* numpy.sum((y-a*x-b)**2); sd=numpy.sqrt(sd)  
   
    return sd  

def predband(xd,yd,a,b,conf=0.95,x=None):
    """
Calculates the prediction band of the linear regression model at the desired confidence
level.
Clarification of the difference between confidence and prediction bands:
"The 2sigma confidence interval is 95% sure to contain the best-fit regression line.
This is not the same as saying it will contain 95% of the data points. The prediction bands are
further from the best-fit line than the confidence bands, a lot further if you have many data
points. The 95% prediction interval is the area in which you expect 95% of all data points to fall."
(from http://graphpad.com/curvefit/linear_regression.htm)
Arguments:
- conf: desired confidence level, by default 0.95 (2 sigma)
- xd,yd: data arrays
- a,b: linear fit parameters as in y=ax+b
- x: (optional) array with x values to calculate the confidence band. If none is provided, will
 by default generate 100 points in the original x-range of the data.
 
Usage:
>>> lpb,upb,x=nemmen.predband(all.kp,all.lg,a,b,conf=0.95)
calculates the prediction bands for the given input arrays
>>> pylab.fill_between(x, lpb, upb, alpha=0.3, facecolor='gray')
plots a shaded area containing the prediction band  
Returns:
Sequence (lpb,upb,x) with the arrays holding the lower and upper confidence bands
corresponding to the [input] x array.
References:
1. http://www.JerryDallal.com/LHSP/slr.htm, Introduction to Simple Linear Regression, Gerard
E. Dallal, Ph.D.
Rodrigo Nemmen
v1 Dec. 2011
v2 Jun. 2012: corrected bug in dy.
    """
    alpha=1.-conf   # significance
    n=xd.size   # data sample size
    if x==None: x=numpy.linspace(xd.min(),xd.max(),100)
    # Predicted values (best-fit model)
    y=a*x+b
    # Auxiliary definitions
    sd=scatterfit(xd,yd,a,b)    # Scatter of data about the model
    sxd=numpy.sum((xd-xd.mean())**2)
    sx=(x-xd.mean())**2 # array
    # Quantile of Student's t distribution for p=1-alpha/2
    q=scipy.stats.t.ppf(1.-alpha/2.,n-2)
    # Prediction band
    dy=q*sd*numpy.sqrt( 1.+1./n + sx/sxd )
    upb=y+dy    # Upper prediction band
    lpb=y-dy    # Lower prediction band
    return lpb,upb,x
