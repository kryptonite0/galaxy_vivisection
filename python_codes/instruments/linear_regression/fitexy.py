import numpy as np

precisionlist = [1., 0.1, 0.01, 0.001, 0.0001]
#aprecision = 0.01
#bprecision = 0.01
epsilonprecision = 0.01

def modfitexy(x,sigx,y,sigy,amin,amax,bmin,bmax,epsilon0):

	N = len(x) 
	#epsilon = epsilon0
	chisqmin = N
	afit,bfit = -9999.,-9999.
	
	for precision in precisionlist:
		
		epsilon = epsilon0		
		while chisqmin>(N-2):  
		
			for a in np.arange(amin,amax,precision):
				
				for b in np.arange(bmin,bmax,precision):
					
					chisq = 0.0
					for i in range(0,N):
						numer = y[i] - (a + b*x[i])
						denom = (sigy[i])**2 + (b**2)*(sigx[i])**2 + epsilon**2
						chisq = chisq + numer**2/denom
					if chisq <=chisqmin:
						chisqmin = chisq
						afit = a
						bfit = b
			
			print 'chisq, chisqmin =', chisq, chisqmin
			epsilon = epsilon + epsilonprecision
			print 'epsilon =', epsilon
		
		amin,amax = (afit-2*precision),(afit+2*precision)	
		bmin,bmax = (bfit-2*precision),(bfit+2*precision)	
					
	#print chisqmin,afit,bfit
	absscat = 0.0
	for i in range(0,N):
		absscat = absscat + (y[i] - afit - bfit*x[i])**2
	absscat = (absscat/N)**0.5 # check if N or N-1	
	
	print '--------------------------'
	print 'y = a + bx'
	print 'Input values: amin =', amin, ' amax =', amax, ' bmin =', bmin, ' bmax =', bmax 
	print 'a =', afit, '; b =', bfit, '; epsilon =', epsilon
	print 'absolute scatter Delta =', absscat					

