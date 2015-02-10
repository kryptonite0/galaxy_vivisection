import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pyfits

from scipy.misc import imread
import matplotlib.cbook as cbook

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 32})

path = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/images/'

def mbhmagkPlot():
	
	fig, ax = plt.subplots()
	ax.set_yscale('log')
        ax.axis([-18.01,-27.99,10**5.5,10**10.3])
        plt.xlabel(r'$M_{\rm K} \rm~[mag]$', labelpad=20)
        plt.ylabel(r'$M_{\rm BH} \rm~[M_\odot]$', labelpad=20)
	plt.subplots_adjust(left=0.20,bottom=0.20)
	
	xxx = np.arange(-29,-17,0.01)
	yyy_c = 10**(9.05 + (-0.44)*(xxx+25))
	yyy_c1 = 10**((9.05+0.09) + (-0.44+0.08)*(xxx+25))
	yyy_c2 = 10**((9.05+0.09) + (-0.44-0.08)*(xxx+25))
	yyy_c3 = 10**((9.05-0.09) + (-0.44+0.08)*(xxx+25))
	yyy_c4 = 10**((9.05-0.09) + (-0.44-0.08)*(xxx+25))
	yyy_cup = yyy_c1
	for i in range(0,len(yyy_cup)):
		if yyy_cup[i] < yyy_c2[i]:
			yyy_cup[i] = yyy_c2[i]
	yyy_clo = yyy_c3
	for i in range(0,len(yyy_clo)):
		if yyy_clo[i] > yyy_c4[i]:
			yyy_clo[i] = yyy_c4[i]
			
	yyy_s = 10**(7.39 + (-1.09)*(xxx+22.5))
	yyy_s1 = 10**((7.39+0.14) + (-1.09+0.22)*(xxx+22.5))
	yyy_s2 = 10**((7.39+0.14) + (-1.09-0.22)*(xxx+22.5))
	yyy_s3 = 10**((7.39-0.14) + (-1.09+0.22)*(xxx+22.5))
	yyy_s4 = 10**((7.39-0.14) + (-1.09-0.22)*(xxx+22.5))
	yyy_sup = yyy_s1
	for i in range(0,len(yyy_sup)):
		if yyy_sup[i] < yyy_s2[i]:
			yyy_sup[i] = yyy_s2[i]
	yyy_slo = yyy_s3
	for i in range(0,len(yyy_slo)):
		if yyy_slo[i] > yyy_s4[i]:
			yyy_slo[i] = yyy_s4[i]

	plt.plot(xxx, yyy_c, color='red', linewidth=3)
	#plt.plot(xxx, yyy_cup, color='red', linewidth=3, ls='--')
	#plt.plot(xxx, yyy_clo, color='red', linewidth=3, ls='--')
	plt.plot(xxx, yyy_s, color='blue', linewidth=3)
	#plt.plot(xxx, yyy_sup, color='blue', linewidth=3, ls='--')
	#plt.plot(xxx, yyy_slo, color='blue', linewidth=3, ls='--')
	
	plt.fill_between(xxx, yyy_cup, yyy_clo, facecolor='red', alpha=0.3)
	plt.fill_between(xxx, yyy_sup, yyy_slo, facecolor='blue', alpha=0.3)
	
	x_n1332 = np.asarray([-24.32])
	y_n1332 = np.asarray([14*10**8])
	yerr_n1332 = np.asarray([2*10**8])
	plt.scatter(x_n1332, y_n1332, c='black', s=80)
        #plt.errorbar(x_n1332, y_n1332, xerr=0, yerr=yerr_n1332, ecolor='black', fmt=None, elinewidth=1.5, capthick=1.5) 
	
	#plt.show()
	plt.savefig(path+'mbhmagk_n1332.pdf', format='pdf', dpi=1000)
	
def main():
	mbhmagkPlot()
	
main()	
