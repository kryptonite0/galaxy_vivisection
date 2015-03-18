import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy
from instruments import b_n
from instruments.linear_regression import colorline
from instruments.linear_regression import akaike



fig, ax = plt.subplots()

#x = np.asarray([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
#y = np.asarray([1.1,1.9,3.1,3.9,5.05,5.99,7.1,7.8,9,10.05,10.55,11,11.45,12.05,12.45,12.9,13.4,13.95])
#
#err_y = y*[0.0] + 0.05
#
#x1 = np.asarray([1,2,3,4,5,6,7,8,9,10])
#y1 = np.asarray([1.1,1.9,3.1,3.9,5.05,5.99,7.1,7.8,9,10.05])
#x2 = np.asarray([11,12,13,14,15,16,17,18])
#y2 = np.asarray([10.55,11,11.45,12.05,12.45,12.9,13.4,13.95])
#
#err_y1 = y1*[0.0] + 0.05
#err_y2 = y2*[0.0] + 0.05

x = np.asarray([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
y = np.asarray([1.1,1.9,3.1,3.9,5.05,5.99,7.1,7.8,9,10.05,11,11.95,13.05,14.1,14.9,16,17,17.95])

err_y = y*[0.0] + 0.05

x1 = np.asarray([1,2,3,4,5,6,7,8,9,10])
y1 = np.asarray([1.1,1.9,3.1,3.9,5.05,5.99,7.1,7.8,9,10.05])
x2 = np.asarray([11,12,13,14,15,16,17,18])
y2 = np.asarray([11,11.95,13.05,14.1,14.9,16,17,17.95])

err_y1 = y1*[0.0] + 0.05
err_y2 = y2*[0.0] + 0.05

print 'all'
A,B,Aerr,Berr,covAB=bces.bces(x,err_y,y,err_y,x*[0.0])
print '---------------------------------'
print 'y = A*(x-<x>) + B '
#print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
#print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
#print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
print '---------------------------------'

m = A[2]
q = B[2]
ax.plot(x,x*m+q, color='k', ls='-', linewidth=2.)

print 'red'
A,B,Aerr,Berr,covAB=bces.bces(x1,err_y1,y1,err_y1,x1*[0.0])
print '---------------------------------'
print 'y = A*(x-<x>) + B '
#print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
#print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
#print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
print '---------------------------------'

m1 = A[2]
q1 = B[2]
ax.plot(x1,x1*m1+q1, color='r', ls='-', linewidth=2.)

print 'blue'
A,B,Aerr,Berr,covAB=bces.bces(x2,err_y2,y2,err_y2,x2*[0.0])
print '---------------------------------'
print 'y = A*(x-<x>) + B '
#print 'OLS(Y|X)    A =', "{0:.4f}".format(A[0]), '+-', "{0:.4f}".format(Aerr[0]), '   B = ', "{0:.4f}".format(B[0]), '+-', "{0:.4f}".format(Berr[0])
#print 'OLS(X|Y)    A =', "{0:.4f}".format(A[1]), '+-', "{0:.4f}".format(Aerr[1]), '   B = ', "{0:.4f}".format(B[1]), '+-', "{0:.4f}".format(Berr[1])
print 'bisector    A =', "{0:.4f}".format(A[2]), '+-', "{0:.4f}".format(Aerr[2]), '   B = ', "{0:.4f}".format(B[2]), '+-', "{0:.4f}".format(Berr[2])
#print 'orthogonal  A =', "{0:.4f}".format(A[3]), '+-', "{0:.4f}".format(Aerr[3]), '   B = ', "{0:.4f}".format(B[3]), '+-', "{0:.4f}".format(Berr[3])
print '---------------------------------'

m2 = A[2]
q2 = B[2]
ax.plot(x2,x2*m2+q2, color='b', ls='-', linewidth=2.)

AICc_single = akaike.get_AICc_singledataset(x,y,err_y,m,q,2)
AICc_double = akaike.get_AICc_doubledataset(x1,y1,err_y1,m1,q1,x2,y2,err_y2,m2,q2,4)

print 
print '-------------------------------'
print 'Akaike result:'
print 'Single power law: AICc =', AICc_single
print 'Double power law: AICc =', AICc_double
print '-------------------------------'

ax.scatter(x1,y1,c='r')
ax.scatter(x2,y2,c='b')
plt.show()
