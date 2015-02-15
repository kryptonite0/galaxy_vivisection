import numpy as np
import bces
import matplotlib.pyplot as plt

tolerance = 0.1

def modfitexy(x,sigx,y,sigy,mode):

	N = len(x)
	# compute first guess of slope and intercept with OLS Y|X method ignoring errors
	cov = [0.0]*x
	# y = A + Bx
	B,A,Berr,Aerr,covAB=bces.bces(x,sigx*[0.0],y,sigy*[0.0],cov)
	index = {'y|x' : 0 , 'x|y': 1}
	
	ainit = A[index[mode]]
	binit = B[index[mode]]
	ainit_err = Aerr[index[mode]]
	binit_err = Berr[index[mode]]

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
	
	
	a, b, epsilon = get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N,mode)
	
	while paramsprecision > tolerance:
		paramsprecision = paramsprecision/3
		a, b, epsilon = get_ab(x,sigx,y,sigy,a,b,paramsprecision,N,mode)
		
	afit, bfit = a, b 	
		
	# print chisqmin,afit,bfit
	# compute absolute scatter
	absscat = 0.0
	for i in range(0,N):
		absscat = absscat + (y[i] - afit - bfit*x[i])**2
	absscat = (absscat/N)**0.5 # check if N or N-1	
	
	print 'y = a + bx'
	print 'Mode fit:', mode
	print 'a =', afit, '; b =', bfit, '; epsilon =', epsilon
	print 'absolute scatter Delta =', absscat
	
	return afit, bfit, epsilon, absscat					

	
def get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N,mode):
	
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
		print 'Epsilon =', epsilon
		print 'Precision =', paramsprecision
		print 
			
		for a in np.arange(ainit-paramsprecision*6,ainit+paramsprecision*6,paramsprecision):
	
			for b in np.arange(binit-paramsprecision*6,binit+paramsprecision*6,paramsprecision):
		
				chisq = 0.000
				for i in range(0,N):
					numer = y[i] - (a + b*x[i])
					if mode == 'y|x':
						denom = (sigy[i])**2 + (b**2)*(sigx[i])**2 + epsilon**2
					elif mode == 'x|y':
						denom = (sigy[i])**2 + (b**2)*(sigx[i])**2 + (b**2)*(epsilon**2)	
					chisq = chisq + numer**2/denom
				if chisq <=chisqmin:
					chisqmin = chisq
					afit = a
					bfit = b
					
	return afit,bfit,epsilon

 
def bisect_modfitexy(x,sigx,y,sigy):
	
	
	afit_1, bfit_1, epsilon_1, absscat_1 = modfitexy(x,sigx,y,sigy,'y|x')
	afit_2, bfit_2, epsilon_2, absscat_2 = modfitexy(x,sigx,y,sigy,'x|y')

	b_bisec = np.tan(0.5*(np.arctan(bfit_1)+np.arctan(bfit_2)))
	a_bisec = (afit_2-afit_1)*(bfit_1-b_bisec)/(bfit_1-bfit_2) + afit_1
	
	print '--------------------------'
	print 'Bisector FITEXY'
	print 'y = a + bx'
	print 'a =', a_bisec
	print 'b =', b_bisec
	
	
        # plot data and regressions
       #fig, ax = plt.subplots()
       #ax.errorbar(x, y, xerr=sigx, yerr=sigy, ecolor='gray', 
       #	fmt='ko', markersize=12, elinewidth=1.2, capthick=1.2, barsabove=False) 
       ##xxx = np.arange(1.33*min(x)-0.33*max(x),1.33*max(x)-0.33*min(x),(max(x)-min(x))/10.)
       #xxx = np.arange(min(x),max(x),(max(x)-min(x))/10.)
       #y_1 = afit_1 + bfit_1*xxx
       #y_2 = afit_2 + bfit_2*xxx
       #y_test = afit_1 - xxx
       #y_bisec = a_bisec + b_bisec*xxx
       #ax.plot(xxx, y_1, ls='--', color='gray', linewidth=2.)
       #ax.plot(xxx, y_2, ls='--', color='gray', linewidth=2.)
       #ax.plot(xxx, y_bisec, ls='-', color='black', linewidth=2.)
       #plt.show()
       #print 'test', [afit_1, afit_2, a_bisec], [bfit_1, bfit_2, b_bisec]
	
	return [afit_1, afit_2, a_bisec], [bfit_1, bfit_2, b_bisec]


