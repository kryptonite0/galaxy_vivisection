import numpy as np

def get_likelihood_singledataset(x,y,err_y,m,q):
	
	L1 = 1
	
	for err_y_i in err_y:
		
		L1 = L1 * (1/(2*np.pi*err_y_i**2))**0.5
	
	lnL1 = np.log(L1)
	
	L2 = 0
	
	for x_i,y_i,err_y_i in zip(x,y,err_y):
	
		model_i = m*x_i + q
		L2 = L2 + 0.5 * ((y_i - model_i)**2/err_y_i**2)
		
	lnL = lnL1 - L2	
	
	return lnL
	
def get_likelihood_doubledataset(x1,y1,err_y1,m1,q1,x2,y2,err_y2,m2,q2):
	
	L1 = 1
	
	for err_y1_i in err_y1:
		L1 = L1 * (1/(2*np.pi*err_y1_i**2))**0.5

	for err_y2_i in err_y2:
		L1 = L1 * (1/(2*np.pi*err_y2_i**2))**0.5
	
	lnL1 = np.log(L1)
	
	L2 = 0
	
	for x1_i,y1_i,err_y1_i in zip(x1,y1,err_y1):
	
		model1_i = m1*x1_i + q1
		L2 = L2 + 0.5 * ((y1_i - model1_i)**2/err_y1_i**2)
		
	for x2_i,y2_i,err_y2_i in zip(x2,y2,err_y2):
	
		model2_i = m2*x2_i + q2
		L2 = L2 + 0.5 * ((y2_i - model2_i)**2/err_y2_i**2)
		
	lnL = lnL1 - L2	
	
	return lnL
	
def get_AICc_singledataset(x,y,err_y,m,q,k):
	
	# m = slope of fit
	# q = intercept of fit
	# k = number of parameters in the model
	
	lnL = get_likelihood_singledataset(x,y,err_y,m,q)
	
	AIC = 2*k - 2*lnL
	
	AICc = AIC + (2*k*(k+1))/(len(x)-k-1)
	
	return AICc
		
def get_AICc_doubledataset(x1,y1,err_y1,m1,q1,x2,y2,err_y2,m2,q2,k):
	
	# m = slope of fit
	# q = intercept of fit
	# k = number of parameters in the model
	
	lnL = get_likelihood_doubledataset(x1,y1,err_y1,m1,q1,x2,y2,err_y2,m2,q2)
	
	AIC = 2*k - 2*lnL
	
	AICc = AIC + (2*k*(k+1))/(len(x1)+len(x2)-k-1)
	
	return AICc
		
