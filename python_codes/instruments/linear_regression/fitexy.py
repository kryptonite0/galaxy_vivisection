import numpy as np
import bces
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import absolutescatter

logfileName = './fitexy.log'

#tolerance = 0.005
tolerance = 0.01  # good one
#tolerance = 0.1

def modfitexy(x,sigx,y,sigy,mode):

	N = len(x)
	# compute first guess of slope and intercept with OLS Y|X method (NON) ignoring errors
	cov = [0.0]*x
	# y = A + Bx
	#B,A,Berr,Aerr,covAB=bces.bces(x,sigx*[0.0],y,sigy*[0.0],cov)
	B,A,Berr,Aerr,covAB=bces.bces(x,sigx,y,sigy,cov)
	index = {'y|x' : 0 , 'x|y': 1}
	
	ainit = A[index[mode]]
	binit = B[index[mode]]
	ainit_err = Aerr[index[mode]]
	binit_err = Berr[index[mode]]

	aprecision = ainit_err/2
	bprecision = binit_err/2
	paramsprecision = min(aprecision,bprecision)
	print '--------------------------'
       #print 'y = a + bx'
       #print 'Inital guess' 
       #print 'a =', ainit, '+-', ainit_err
       #print 'b =', binit, '+-', binit_err
       #print 'precision = ', paramsprecision
       #print '--------------------------'
       #print
	
	
	a, b, epsilon, chisqmin = get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N,mode)
	
	while paramsprecision > tolerance:
		paramsprecision = paramsprecision/3
		a, b, epsilon, chisqmin = get_ab(x,sigx,y,sigy,a,b,paramsprecision,N,mode)
		
	afit, bfit, epsilonfit, chisqminfit = a, b, epsilon, chisqmin 	
	
	epsilonfit_minus, epsilonfit_plus = get_err_epsilon(x,sigx,y,sigy,afit,bfit,epsilonfit,chisqminfit,0.001,N,mode)
	
	perr_epsilonfit = epsilonfit_plus - epsilonfit
	merr_epsilonfit = epsilonfit - epsilonfit_minus
	
	print 'get_ab completed'				
	
	merr_afit,perr_afit,merr_bfit,perr_bfit = get_errors_ab(x,sigx,y,sigy,afit,bfit,ainit_err,binit_err,epsilonfit,chisqmin,tolerance,N,mode)
	print 'get_errors_ab completed'				
		
	# print chisqmin,afit,bfit
	# compute absolute scatter
	absscat = absolutescatter.get_absscatter(x, y, afit, bfit)
	
	if mode == 'x|y':
		epsilonfit = epsilonfit*abs(bfit)
		perr_epsilonfit = perr_epsilonfit*abs(bfit)
		merr_epsilonfit = merr_epsilonfit*abs(bfit)
	
	print 'y = a + bx'
	print 'Mode fit:', mode
	print 'a =', "{0:.2f}".format(afit), '^{+', "{0:.2f}".format(perr_afit), '}_{-', "{0:.2f}".format(merr_afit), '}'
	print 'b =', "{0:.2f}".format(bfit), '^{+', "{0:.2f}".format(perr_bfit), '}_{-', "{0:.2f}".format(merr_bfit), '}' 
	print 'epsilon =', "{0:.2f}".format(epsilonfit), '^{+', "{0:.2f}".format(perr_epsilonfit), '}_{-', "{0:.2f}".format(merr_epsilonfit), '}' 
	print 'absolute scatter Delta =', "{0:.2f}".format(absscat)
	
	return afit,merr_afit,perr_afit, bfit,merr_bfit,perr_bfit, epsilonfit,merr_epsilonfit,perr_epsilonfit, absscat					

	
def get_ab(x,sigx,y,sigy,ainit,binit,paramsprecision,N,mode):
	
	chisqmin = N + 0.000
	afit,bfit = -9999.,-9999.
	epsilon = -paramsprecision
	iteration = -1			
	while chisqmin>(N-2.): 
	
		epsilon = epsilon + paramsprecision
		iteration = iteration + 1
		#print 'Iteration', iteration, '; chisq =', chisqmin, '; N-2 =', (N-2)
	        #print 'Range for a:', ainit-paramsprecision*6,ainit+paramsprecision*6
	        #print 'Range for b:', binit-paramsprecision*6,binit+paramsprecision*6
	        #print 'Epsilon =', epsilon
	        #print 'Precision =', paramsprecision
	        #print 
			
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
	return afit,bfit,epsilon,chisqmin

def get_err_epsilon(x,sigx,y,sigy,afit,bfit,epsilonfit,chisqminfit,precision,N,mode):
	
	epsilon_minus, epsilon_plus = -9999.,-9999.
	
	chisqmin = chisqminfit
	epsilon = epsilonfit
	
	while chisqmin<((N-2.)*(1.+(2./N)**0.5)): 
	
		epsilon = epsilon - precision
		#print 'e', epsilon
			
		chisq = 0.000
		for i in range(0,N):
			numer = y[i] - (afit + bfit*x[i])
			if mode == 'y|x':
				denom = (sigy[i])**2 + (bfit**2)*(sigx[i])**2 + epsilon**2
			elif mode == 'x|y':
				denom = (sigy[i])**2 + (bfit**2)*(sigx[i])**2 + (bfit**2)*(epsilon**2)	
			chisq = chisq + numer**2/denom
		chisqmin = chisq
		#print 'chi', chisqmin, ((N-2.)*(1.+(2./N)**0.5))
		print 'chisq', chisqmin
		print 'e', epsilon
		epsilon_minus = epsilon
	
	chisqmin = chisqminfit
	epsilon = epsilonfit
	
	while chisqmin>((N-2.)*(1.-(2./N)**0.5)): 
	
		epsilon = epsilon + precision
			
		chisq = 0.000
		for i in range(0,N):
			numer = y[i] - (afit + bfit*x[i])
			if mode == 'y|x':
				denom = (sigy[i])**2 + (bfit**2)*(sigx[i])**2 + epsilon**2
			elif mode == 'x|y':
				denom = (sigy[i])**2 + (bfit**2)*(sigx[i])**2 + (bfit**2)*(epsilon**2)	
			chisq = chisq + numer**2/denom
		chisqmin = chisq
		print 'chisq', chisqmin
		print 'e', epsilon
		epsilon_plus = epsilon
	
	return epsilon_minus, epsilon_plus

 
def bisect_modfitexy(x,sigx,y,sigy):

	logfile = open(logfileName, 'a')
	
	afit_1,merr_afit_1,perr_afit_1, bfit_1,merr_bfit_1,perr_bfit_1, epsilonfit_1,merr_epsilonfit_1,perr_epsilonfit_1, absscat_1 = modfitexy(x,sigx,y,sigy,'y|x')
	logfile.write('--------------------------------------------------------- \n')
	logfile.write('\ny = a + bx \n\nMode fit: y|x \na =' + str(afit_1) + '+' + str(perr_afit_1) + '-' + str(merr_afit_1) + '\n')
	logfile.write('b =' + str(bfit_1) + '+' + str(perr_bfit_1) + '-' + str(merr_bfit_1) + '\n')
	logfile.write('epsilon =' + str(epsilonfit_1) + '+' + str(perr_epsilonfit_1) + '-' + str(merr_epsilonfit_1) + '\n')
	logfile.write('absolute scatter Delta =' + str(absscat_1) + '\n')
	afit_2,merr_afit_2,perr_afit_2, bfit_2,merr_bfit_2,perr_bfit_2, epsilonfit_2,merr_epsilonfit_2,perr_epsilonfit_2, absscat_2 = modfitexy(x,sigx,y,sigy,'x|y')
	logfile.write('\nMode fit: x|y \na =' + str(afit_2) + '+' + str(perr_afit_2) + '-' + str(merr_afit_2) + '\n')
	logfile.write('b =' + str(bfit_2) + '+' + str(perr_bfit_2) + '-' + str(merr_bfit_2) + '\n')
	logfile.write('epsilon =' + str(epsilonfit_2) + '+' + str(perr_epsilonfit_2) + '-' + str(merr_epsilonfit_2) + '\n')
	logfile.write('absolute scatter Delta =' + str(absscat_2))
	
	b_bisec = np.tan(0.5*(np.arctan(bfit_1)+np.arctan(bfit_2)))
	a_bisec = (afit_2-afit_1)*(bfit_1-b_bisec)/(bfit_1-bfit_2) + afit_1
	
	## compute errors
	# error on slope
	perr_alpha_1 = np.arctan(bfit_1+perr_bfit_1) - np.arctan(bfit_1)
	merr_alpha_1 = np.arctan(bfit_1) - np.arctan(bfit_1-merr_bfit_1)
	perr_alpha_2 = np.arctan(bfit_2+perr_bfit_2) - np.arctan(bfit_2)
	merr_alpha_2 = np.arctan(bfit_2) - np.arctan(bfit_2-merr_bfit_2)
	
	perr_alpha_bisec = (perr_alpha_1**2+perr_alpha_2**2)**0.5/2**0.5
	merr_alpha_bisec = (merr_alpha_1**2+merr_alpha_2**2)**0.5/2**0.5
	
	perr_b_bisec = np.tan(perr_alpha_bisec+np.arctan(b_bisec)) - b_bisec
	merr_b_bisec = b_bisec - np.tan(np.arctan(b_bisec)-merr_alpha_bisec)
	
	#error on intercept
	#true only in the approximation afit_1 ~ afit_2
	perr_a_bisec = (perr_afit_1**2 + perr_afit_2**2)**0.5/2**0.5
	merr_a_bisec = (merr_afit_1**2 + merr_afit_2**2)**0.5/2**0.5
	
	absscat_bisec = absolutescatter.get_absscatter(x, y, a_bisec, b_bisec)

	
	print '--------------------------'
	print 'Bisector FITEXY'
	print 'y = a + bx'
	print 'a =', "{0:.2f}".format(a_bisec), '^{+', "{0:.2f}".format(perr_a_bisec), '}_{-', "{0:.2f}".format(merr_a_bisec), '}' 
	print 'b =', "{0:.2f}".format(b_bisec), '^{+', "{0:.2f}".format(perr_b_bisec), '}_{-', "{0:.2f}".format(merr_b_bisec), '}' 
	print 'absolute scatter Delta =', "{0:.2f}".format(absscat_bisec)
	
	logfile.write('\n\nBisector FITEXY \n')
	logfile.write('a =' + str(a_bisec) + '+' + str(perr_a_bisec) + '-' + str(merr_a_bisec) + '\n')
	logfile.write('b =' + str(b_bisec) + '+' + str(perr_b_bisec) + '-' + str(merr_b_bisec) + '\n')
	logfile.write('--------------------------------------------------------- \n')
	logfile.close()
	
	#epsilon_bisec = get_epsilon_bisec(x,sigx,y,sigy,a_bisec,b_bisec,paramsprecision,N)
	#print 'epsilon =', epsilon_bisec
	
	
       ## plot data and regressions
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
       ##ax.plot(xxx, y_bisec, ls='-', color='black', linewidth=2.)
       #plt.show()
       ##print 'test', [afit_1, afit_2, a_bisec], [bfit_1, bfit_2, b_bisec]
	
	#return [afit_1, afit_2, a_bisec], [bfit_1, bfit_2, b_bisec]
	return a_bisec,perr_a_bisec,merr_a_bisec,b_bisec,perr_b_bisec,merr_b_bisec


def get_errors_ab(x,sigx,y,sigy,afit,bfit,ainit_err,binit_err,epsilonfit,chisqmin,tolerance,N,mode):
	
	a_range = np.arange(afit-3*ainit_err,afit+3*ainit_err+tolerance,tolerance)
	b_range = np.arange(bfit-3*binit_err,bfit+3*binit_err+tolerance,tolerance)
	chisqgrid = np.zeros((len(a_range), len(b_range)))
		
	for n1 in range(len(a_range)):
		
		for n2 in range(len(b_range)):
			
			chisq = 0.000
			for i in range(0,N):
				numer = y[i] - (a_range[n1] + b_range[n2]*x[i])
				if mode == 'y|x':
					denom = (sigy[i])**2 + (b_range[n2]**2)*(sigx[i])**2 + epsilonfit**2
				elif mode == 'x|y':
					denom = (sigy[i])**2 + (b_range[n2]**2)*(sigx[i])**2 + (b_range[n2]**2)*(epsilonfit**2)
				chisq = chisq + numer**2/denom
			
			chisqgrid[n1,n2] = chisq
			
	print 'chisq grid calculated'
	
        chisqimg = np.zeros((len(a_range), len(b_range))) + 1
        chisqimg[chisqgrid>chisqmin+1] = 0
        
        chisqimg_collapsed_a = np.sum(chisqimg, axis=1)
        chisqimg_collapsed_b = np.sum(chisqimg, axis=0)
        
       #print 'producing ab plot...'
       #fig, ax = plt.subplots()
       #ax.imshow(chisqimg, interpolation='none', cmap=cm.Greys_r)
       #ax.axhline(y=0.5*len(chisqimg_collapsed_a)-1, linewidth=2, color = 'red')
       #ax.axvline(x=0.5*len(chisqimg_collapsed_b)-1, linewidth=2, color = 'red')
       #plt.show()
	
	for j in range(len(chisqimg_collapsed_a)-1):
		if chisqimg_collapsed_a[j]>0 and chisqimg_collapsed_a[j-1]<0.5:
			jleft = j
		if chisqimg_collapsed_a[j+1]<0.5 and chisqimg_collapsed_a[j]>0:	
			jright = j
			break
			
			
	for k in range(len(chisqimg_collapsed_b)-1):
		if chisqimg_collapsed_b[k]>0 and chisqimg_collapsed_b[k-1]<0.5:
			kleft = k
		if chisqimg_collapsed_b[k+1]<0.5 and chisqimg_collapsed_b[k]>0:	
			kright = k
			break
			
	
	merr_afit = afit - a_range[jleft]
	perr_afit = a_range[jright] - afit	
	merr_bfit = bfit - b_range[kleft]
	perr_bfit = b_range[kright] - bfit
	
	err_afit = (a_range[jright] - a_range[jleft])/2
	err_bfit = (b_range[kright] - b_range[kleft])/2	
			
       #fig, ax = plt.subplots()
       #ax.imshow(chisqimg, interpolation='none', cmap=cm.Greys_r)
       #ax.axhline(y=jleft, linewidth=2, color = 'red')
       #ax.axhline(y=jright, linewidth=2, color = 'red')
       #ax.axhline(y=0.5*len(chisqimg_collapsed_a)-1, linewidth=2, color = 'red')
       #ax.axvline(x=kleft, linewidth=2, color = 'red')
       #ax.axvline(x=kright, linewidth=2, color = 'red')
       #ax.axvline(x=0.5*len(chisqimg_collapsed_b)-1, linewidth=2, color = 'red')
       #
       #plt.show()
	#print 'errors', merr_afit,perr_afit,merr_bfit,perr_bfit
	#print 'errors', err_afit,err_bfit
	return merr_afit,perr_afit,merr_bfit,perr_bfit
	
	






