import numpy as np
import bces

tolerance = 0.001

def modfitexy(x,sigx,y,sigy):

	N = len(x)
	# compute first guess of slope and intercept with OLS Y|X method
	cov = [0.0]*x
	# y = A + Bx
	B,A,Berr,Aerr,covAB=bces.bces(x,sigx,y,sigy,cov)
	ainit = A[0]
	binit = B[0]
	ainit_err = Aerr[0]
	binit_err = Berr[0]

	aprecision = ainit_err/2
	bprecision = binit_err/2
	paramsprecision = min(aprecision,bprecision)
	print '--------------------------'
	print 'y = a + bx'
	print 'Inital guess' 
	print 'a =', ainit, '+-', ainit_err
	print 'b =', binit, '+-', binit_err
	print 'precision = ', paramsprecision
	print '--------------------------'
	print
	
	
	a, b, epsilon = get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N)
	
	while paramsprecision > tolerance:
		paramsprecision = paramsprecision/3
		a, b, epsilon = get_ab(x,sigx,y,sigy,a,b,paramsprecision,N)
		
	afit, bfit = a, b 	
		
	# print chisqmin,afit,bfit
	# compute absolute scatter
	absscat = 0.0
	for i in range(0,N):
		absscat = absscat + (y[i] - afit - bfit*x[i])**2
	absscat = (absscat/N)**0.5 # check if N or N-1	
	
	print 'a =', afit, '; b =', bfit, '; epsilon =', epsilon
	print 'absolute scatter Delta =', absscat					

	
def get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N):
	
	chisqmin = N + 0.000
	afit,bfit = -9999.,-9999.
	epsilon = -paramsprecision
	iteration = -1			
	while chisqmin>(N-2): 
	
		epsilon = epsilon + paramsprecision
		iteration = iteration + 1
		print 'Iteration', iteration, '; chisq =', chisqmin, '; N-2 =', (N-2)
		print 'Range for a:', ainit-paramsprecision*6,ainit+paramsprecision*6
		print 'Range for b:', binit-paramsprecision*6,binit+paramsprecision*6
		print 'Precision =', paramsprecision
		print 
			
		for a in np.arange(ainit-paramsprecision*6,ainit+paramsprecision*6,paramsprecision):
	
			for b in np.arange(binit-paramsprecision*6,binit+paramsprecision*6,paramsprecision):
		
				chisq = 0.000
				for i in range(0,N):
					numer = y[i] - (a + b*x[i])
					denom = (sigy[i])**2 + (b**2)*(sigx[i])**2 + epsilon**2
					chisq = chisq + numer**2/denom
				if chisq <=chisqmin:
					chisqmin = chisq
					afit = a
					bfit = b
					
	return afit,bfit,epsilon
	
def get_ep	
