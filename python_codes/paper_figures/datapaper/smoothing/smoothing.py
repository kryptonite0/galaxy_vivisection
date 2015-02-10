import numpy as np
#from lmfit import minimize, Parameters, report_fit
#import inspect
#from scipy.ndimage.filters import gaussian_filter1d, correlate1d
#from scipy.special import gammainc, gamma
#import pylab
import matplotlib.pyplot as plt
#import shutil
import matplotlib

from ellipseOutput import readEllipseOutput

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 32})


    
def readAndPlot(ellipseOutput):

	#build data arrays

	datatab = readEllipseOutput(ellipseOutput) 

	#exclude R=0 point
	datatab = datatab[datatab['sma']>0]
	#exclude data fainter than (3) sky RMS
	datatab = datatab[datatab['intens']>skyRMS]

	row = datatab['row']
	sma = datatab['sma']
	intens = datatab['intens']
	    
	rrr = sma * 1.22
	mu = 17.26 - 2.5*np.log10(intens) + 2.5*np.log10(1.22**2)

		    
	return rrr, mu
	
####################################################

skyRMS = 3*0.01

#smoothedEllipseOutput = 'm31_op_overmos_sky_gauss60_combscale.ell'
#noisyEllipseOutput = 'm31_op_overmos_sky_linscale1px.ell'   
smoothedEllipseOutput = 'm81_op_overmos_sky_gauss5_linscale1px.ell'
noisyEllipseOutput = 'm81_op_overmos_sky_linscale1px.ell'   
                
    

fig, ax1 = plt.subplots()

rrr_smoo, mu_smoo = readAndPlot(smoothedEllipseOutput)
rrr_nois, mu_nois = readAndPlot(noisyEllipseOutput)
maxrrr = max(rrr_smoo)
minmu = min(mu_smoo)
maxmu = max(mu_smoo)

plt.axis([-100,maxrrr+100,maxmu+1,minmu-1])
#plt.plot(rrr_smoo, mu_smoo, c='blue', linewidth=3)
#plt.plot(rrr_nois, mu_nois, c='red', linewidth=3, ls='--')
plt.scatter(rrr_smoo, mu_smoo, c='blue')
plt.scatter(rrr_nois, mu_nois, c='red')



plt.xlabel(r'$R \rm~[arcsec]$', labelpad=20)
plt.ylabel(r'$\mu$ [mag arcsec$^{-2}$]', labelpad=20)
#plt.title(title)
plt.subplots_adjust(left=0.20,bottom=0.20)
	

plt.show()
#bestfitFig.subplots_adjust(wspace=0, hspace=0)
#bestfitFig.savefig('m31.eps', format='eps', dpi=1000)


########## FIT AND PLOT ################

	



