import numpy as np
from instruments.linear_regression import fitexy

#data = open('test.dat')
#data = open('Savorgnan_2014.dat')
#data = open('msigma.out')
#data = open('allcoreell.dat')
#data = open('goodnocoresp.dat')
data = open('allnocorespi.dat')

xlist = []
sigxlist = []
ylist = []
sigylist = []

for line in data:
	if line[0] != '#':
		xlist.append(float(line.split()[0]))
		sigxlist.append(float(line.split()[1]))	
		ylist.append(float(line.split()[2]))
		sigylist.append(float(line.split()[3]))	
	

x = np.asarray(xlist)
sigx = np.asarray(sigxlist)
y = np.asarray(ylist)
sigy = np.asarray(sigylist)

#fitexy.modfitexy(x,sigx,y,sigy,'y|x')
#fitexy.modfitexy(x,sigx,y,sigy,'x|y')

fitexy.bisect_modfitexy(x,sigx,y,sigy)

