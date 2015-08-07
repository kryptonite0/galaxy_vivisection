import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 26})


path_paper_figures = '/Users/gsavorgnan/galaxy_vivisection/papers/lentiptical/images/'

def spiral(ax, color, x0, y0, maxxx):
	t1 = np.arange(0.5,maxxx,0.01)
	t2 = np.arange(0.5,maxxx,0.01)
	x1 = t1*np.cos(t1)
	y1 = t1*np.sin(t1)
	x2 = -t2*np.cos(t2)
	y2 = -t2*np.sin(t2)
	
	t = np.arange(0,2*np.pi,0.01)
	x = 1.1*np.cos(t)
	y = 1*np.sin(t)
	
	ax.plot(x0+x1, y0+y1, color=color, lw=6)
	ax.plot(x0+x2, y0+y2, color=color, lw=6)
	#ax.plot(x0+x,  y0+y,  color=color, lw=2)
	#ax.fill_between((x0+x), (y0+0), (y0+y), edgecolor='', facecolor=color)







def Sa():
	
        fig, ax = plt.subplots()
	
	spiral(ax, 'blue', 1, 1, 14)

        plt.axis([-30,30,-30,30])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
        #plt.show()
	plt.savefig(path_paper_figures + 'Sa.pdf', format='pdf', dpi=1000)

def Sb():
	
        fig, ax = plt.subplots()
	
	spiral(ax, 'blue', 1, 1, 8)

        plt.axis([-20,20,-20,20])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
        #plt.show()
	plt.savefig(path_paper_figures + 'Sb.pdf', format='pdf', dpi=1000)

def Sc():
	
        fig, ax = plt.subplots()
	
	spiral(ax, 'blue', 1, 1, 6)

        plt.axis([-15,15,-15,15])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
        #plt.show()
	plt.savefig(path_paper_figures + 'Sc.pdf', format='pdf', dpi=1000)

def Sd():
	
        fig, ax = plt.subplots()
	
	spiral(ax, 'blue', 1, 1, 5)

        plt.axis([-10,10,-10,10])
        #plt.xlabel(r'$M_{\rm *,sph}\rm~[M_\odot]$', labelpad=13)
        #plt.ylabel(r'$M_{\rm BH} \rm ~[M_\odot]$', labelpad=13)
        #plt.show()
	plt.savefig(path_paper_figures + 'Sd.pdf', format='pdf', dpi=1000)



Sa()
Sb()
Sc()
Sd()
