import numpy as np
import fitexy

#data = open('test.dat')
data = open('Savorgnan_2014.dat')

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

fitexy.modfitexy(x,sigx,y,sigy,2.,8.,0.,4.,0.0)
#fitexy.modfitexy(x,sigx,y,sigy,5.,9.,1,10.,0.0)
