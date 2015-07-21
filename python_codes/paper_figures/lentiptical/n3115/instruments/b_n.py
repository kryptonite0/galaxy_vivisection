import numpy as np
from scipy.special import gammainc
import math
from decimal import Decimal
from scipy.optimize import minimize_scalar

######## SETTINGS ##########################

precisionResid = 0.000001

######## COMPUTE b_n #######################

def b_nResid(k, x):
        
        out = gammainc(k, x) - 0.50000000000000
        return out

def computeb_n(n):

        k = 2*n
		
	b=2*n - 0.333
	if (b<0):
		b=0.000
		
	zero  = b_nResid(k, b)	
	increment = 1.000
	
	while (abs(zero) > precisionResid):
		
		zero  = b_nResid(k, b)
		
		bplus = b + increment
		bminus = b - increment
		if (bminus<0):
			bminus = 0.000
		
		zeroplus = b_nResid(k, bplus)	
		zerominus = b_nResid(k, bminus)
		
		if (abs(zeroplus) < precisionResid):
			b = bplus
			break
       		elif (abs(zerominus) < precisionResid):
			b = bminus
			break
       		else:
			if ( abs(zeroplus) < abs(zero) ):
				b = bplus
			elif ( abs(zerominus) < abs(zero) ):
				b = bminus
			increment = increment/2.	
				
        return b
	




