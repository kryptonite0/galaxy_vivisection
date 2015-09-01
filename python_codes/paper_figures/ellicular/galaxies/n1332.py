# at the moment, works well in pixels

import numpy as np
from lmfit import minimize, Parameters, report_fit
import inspect
from scipy.ndimage.filters import gaussian_filter1d, correlate1d
from scipy.special import gammainc, gamma
import sys
import pylab
import matplotlib.pyplot as plt
import shutil

# local folders of project

from reading.dataFileNames import readDataFileNames
from reading.ellipseOutput import readEllipseOutput
from reading.inputModel import readInputModel
from reading.excludeData import readExcludedData
from reading.bestFitModel import readBestFitModel

from cls.cls import *
from cls.settings import *
from cls.psf import *

from instruments.fitResult import printFitResult
from instruments.model import buildModel
from instruments.output import *

from instruments.residual import residual
from instruments import b_n
#from instruments.plotting import *

from modelComponents.disc import buildDisc
from modelComponents.ferrer import buildFerrer
from modelComponents.gaussian import buildGaussian
from modelComponents.psf import buildPsf
from modelComponents.psfWing import buildPsfWing
from modelComponents.sersic import buildSersic
from modelComponents.truncDisc import buildTruncDisc
from modelComponents.truncSersic import buildTruncSersic
from modelComponents.gaussRing import buildGaussRing
from instruments.model import buildModel

from convolution.psfConvolution import doPsfConvolution


from matplotlib.pyplot import fill_between
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc, font_manager

from textwrap import wrap



#################################################### 

moffatPsf = PsfFunction()
moffatPsf.name = 'moffat'
moffatPsf.moffatAlpha = (1.61467/(2*np.sqrt(2**(1/4.39)-1)) ) * Settings.pxlToArcsec  # alpha = fwhm / (2 * sqrt(2**(1/beta) - 1) )
moffatPsf.moffatBeta = 4.39  
        
psfList = [moffatPsf]
	
gaussianSmoothing = PsfFunction()
gaussianSmoothing.name = 'gaussian'
gaussianSmoothing.gaussianFWHM = 5 * 2.3548 * Settings.pxlToArcsec 
        
####################################################

terminal = sys.stdout


########## FIT AND PLOT ################

def readFitAndPlot(excludedRangeList, galaxy, axisFit, psfFunction, sampling, bestfitFig, equivalentAxisFit):

    #build data arrays
    
    datatab = readEllipseOutput(ellipseOutput) 
    
    #exclude R=0 point
    datatab = datatab[datatab['sma']>0]
    #exclude data fainter than (3) sky RMS
    datatab = datatab[datatab['intens']>skyRMS]
    
    row = datatab['row']
    sma = datatab['sma']
    intens = datatab['intens']
    intens_err = datatab['intens_err']
    #pix_var = datatab['pix_var']
    #rms = datatab['rms']
    ellip = datatab['ellip']
    ellip_err = datatab['ellip_err']
    #PA = datatab['PA']
    #PA_err = datatab['PA_err']
    
    #ellip[ellip=='nan'] = 0 # !!! check this
        
    rrr = sma * Settings.pxlToArcsec
    mu = Settings.zeropoint - 2.5*np.log10(intens) + 2.5*np.log10(Settings.pxlToArcsec**2)
    
    maxsma_arcsec = max(sma) * Settings.pxlToArcsec
    minmu = min(mu)
    maxmu = max(mu)
    
    # build parameters and component list (not really the best way to do it)
    
    componentslist, params, psfwing_02pxscale_datatab, psfwing_logscale_datatab = readInputModel(inputModel, equivalentAxisFit, Settings)
    
    # exclude data
    
    goodIndexes = list(range(0,len(rrr),1))
    badIndexes = list()
    
    if excludedRangeList is not None:
        for r in rrr:
            index = list(rrr).index(r)
            for excludedRange in excludedRangeList:
                x1 = min(excludedRange)
                x2 = max(excludedRange)
                if (r>x1 and r<x2):
                    if index in goodIndexes:
                        goodIndexes.remove(index)
                    if index not in badIndexes:
                        badIndexes.append(index)
    
    goodIndexes_integer = goodIndexes
    badIndexes_integer = badIndexes		

    if (sampling != 'log'):
        if not Settings.smoothing:
            for s in sma:
                if not (s == int(s)):
                    index = list(sma).index(s)
                    if index in goodIndexes_integer:
                        goodIndexes_integer.remove(index)
                    elif index in badIndexes_integer:
                        badIndexes_integer.remove(index)	
        elif Settings.smoothing:
            for s in sma:
                if not (s%(int(gaussianSmoothing.gaussianFWHM / (2.35 * Settings.pxlToArcsec))) == 0):
                    index = list(sma).index(s)
                    if index in goodIndexes_integer:
                        goodIndexes_integer.remove(index)
                    elif index in badIndexes_integer:
                        badIndexes_integer.remove(index)	
                            

    good_rrr = np.array([])
    bad_rrr = np.array([])
    good_intens = np.array([])
    bad_intens = np.array([])
    good_intens_err = np.array([])
    bad_intens_err = np.array([])
    good_ellip = np.array([])
    bad_ellip = np.array([])
    
    for i in range(0,len(rrr),1):
        if i in goodIndexes:
            good_rrr = np.append(good_rrr, rrr[i])
            good_intens = np.append(good_intens, intens[i])
            good_intens_err = np.append(good_intens_err, intens_err[i])
            good_ellip = np.append(good_ellip, ellip[i])
        elif i in badIndexes:
            bad_rrr = np.append(bad_rrr, rrr[i])
            bad_intens = np.append(bad_intens, intens[i])
            bad_intens_err = np.append(bad_intens_err, intens_err[i])	
            bad_ellip = np.append(bad_ellip, ellip[i])
    
    good_mu = Settings.zeropoint - 2.5*np.log10(good_intens) + 2.5*np.log10(Settings.pxlToArcsec**2)
    bad_mu = Settings.zeropoint - 2.5*np.log10(bad_intens) + 2.5*np.log10(Settings.pxlToArcsec**2)  
    good_mu_up_err = 2.5*np.log10(1 - good_intens_err/good_intens)
    good_mu_lo_err = 2.5*np.log10(1 + good_intens_err/good_intens)
    bad_mu_up_err = 2.5*np.log10(1 - bad_intens_err/bad_intens)
    bad_mu_lo_err = 2.5*np.log10(1 + bad_intens_err/bad_intens)
    
        
    if (equivalentAxisFit):
        good_rrr = good_rrr * np.sqrt(1. - good_ellip) # semi equivalent axis	
        bad_rrr = bad_rrr * np.sqrt(1. - bad_ellip) 
        rrr = rrr * np.sqrt(1. - ellip) 
    
    
    #bestFitParamsFileName = galaxy + '_' + axisFit + '_' + gaussianPsf.name + '_comb_par_SM.dat'
    bestFitParamsFileName = galaxy + '_' + axisFit + '_' + moffatPsf.name + '_comb_par_SM.dat'
    
    componentslist, finalparams = readBestFitModel(bestFitParamsFileName)
    
    bestFitParamsFile = open(bestFitParamsFileName)
    lines = bestFitParamsFile.readlines()
    for line in lines:
        if (line.split()[0] == 'Delta'):                
            deltarms = float(line.split()[1])
    bestFitParamsFile.close()        

    y_finalmodel_sb, good_y_finalmodel_sb = buildModel(finalparams, componentslist, rrr, goodIndexes, psfFunction, Settings.plotConvolvedFinalModel, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing)
    resid = mu - y_finalmodel_sb
    good_resid = good_mu - good_y_finalmodel_sb


#### produce figures here	

    bestfitFig = createFigure(bestfitFig, psfFunction, equivalentAxisFit, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, resid, good_resid, componentslist, finalparams, Settings, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing, datatab)



def doFitPanel(fitPanel, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, componentslist, finalparams, Settings, psfFunction, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing):
	
	sensitivity = Settings.zeropoint - 2.5*np.log10(skyRMS) + 2.5*np.log10(Settings.pxlToArcsec**2)
	
	fitPanel.set_xlim([maxsma_arcsec*(-0.1),1.1*maxsma_arcsec])
	fitPanel.set_ylim([maxmu+1,minmu-1])

	if (sensitivity>maxmu+0.5):
		fitPanel.set_ylim([sensitivity+0.5,minmu-1])
	
	minorLocator   = MultipleLocator(xaxisMinorLocator)
	fitPanel.xaxis.set_minor_locator(minorLocator)
#  	majorLocator   = MultipleLocator(xaxisMajorLocator)
#  	fitPanel.xaxis.set_major_locator(majorLocator)

	majorLocator   = MultipleLocator(2)
	minorLocator   = MultipleLocator(1)
	fitPanel.yaxis.set_major_locator(majorLocator)
	fitPanel.yaxis.set_minor_locator(minorLocator)
	
	start, end = fitPanel.get_ylim()
	ticks = np.arange(int(start)-1, int(end), -2)
	fitPanel.yaxis.set_ticks(ticks)
	fitPanel.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
	
	#xxx = np.arange(min(rrr),max(rrr),0.1)
	xxx = rrr
	
	y_finalmodel_sb_plot, good_y_finalmodel_sb_plot = buildModel(finalparams, componentslist, xxx, goodIndexes, psfFunction, Settings.plotConvolvedFinalModel, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing)
	
	fitPanel.plot(xxx, y_finalmodel_sb_plot, 'black')

	for component in componentslist:
		if (component.name == 'sersic'):
			y_component = buildSersic(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)
		elif (component.name == 'disc'):
			y_component = buildDisc(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)			
		elif (component.name == 'gaussian'):
			y_component = buildGaussian(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)	
		elif (component.name == 'psf'):
			y_component = buildPsf(xxx, finalparams, component, psfFunction, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)			
		elif (component.name == 'ferrer'):
			y_component = buildFerrer(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)	
		elif (component.name == 'tdisc'):
			y_component = buildTruncDisc(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)	
		elif (component.name == 'tsersic'):
			y_component = buildTruncSersic(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)	
		elif (component.name == 'gring'):
			y_component = buildGaussRing(xxx, finalparams, component, psfFunction, Settings.plotConvolvedComponent, Settings.zeropoint)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)	
		elif (component.name == 'psfwing'):
			y_component = buildPsfWing(xxx, finalparams, component, sampling, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings)   
# 			y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)			
		else:
			y_component = None
			print 'Not valid function name: ', component.name
		
		if Settings.smoothing:
			y_component = doPsfConvolution(y_component, gaussianSmoothing, xxx)
		y_component_sb = Settings.zeropoint - 2.5*np.log10(y_component)
		fitPanel.plot(xxx, y_component_sb, color=Settings.colordict[component.name])	
	
	if (sampling == 'log'):
		fitPanel.plot(good_rrr, good_mu, 'ko', markersize=4)
		fitPanel.plot(bad_rrr, bad_mu, 'wo', markersize=4)
	elif (sampling == 'comb'):
		good_rrr_undersampled = np.array([])
		good_mu_undersampled = np.array([])
		bad_rrr_undersampled = np.array([])
		bad_mu_undersampled = np.array([])
		for i in range (0,len(rrr),1):
			if i in goodIndexes_integer:
				good_rrr_undersampled = np.append(good_rrr_undersampled, rrr[i])
				good_mu_undersampled = np.append(good_mu_undersampled, mu[i])
			if i in badIndexes_integer:
				bad_rrr_undersampled = np.append(bad_rrr_undersampled, rrr[i])
				bad_mu_undersampled = np.append(bad_mu_undersampled, mu[i])
		fitPanel.plot(good_rrr_undersampled, good_mu_undersampled, 'ko', markersize=4)
		fitPanel.plot(bad_rrr_undersampled, bad_mu_undersampled, 'wo', markersize=4)
		
	fitPanel.axhline(sensitivity, linewidth=2, color='grey', ls = 'dashed')
	
	fitPanel.text(0, (sensitivity-0.3), r'$\Delta = $' + ("{0:.4f}".format(deltarms)))	
	
	if (componentslist[0].number == 1 and componentslist[0].name == 'sersic'):
		r_e = finalparams['par1_1'].value
		mu_e = finalparams['par2_1'].value
		n = finalparams['par3_1'].value		
		#fitPanel.text(maxsma_arcsec*0.7, (minmu+0.5), r'$R_{\rm e} = $' + ("{0:.1f}".format(r_e)), color='red')
		#fitPanel.text(maxsma_arcsec*0.7, (minmu+1.5), r'$\mu_{\rm e} = $' + ("{0:.2f}".format(mu_e)), color='red')
		#fitPanel.text(maxsma_arcsec*0.7, (minmu+2.5), r'$n = $' + ("{0:.1f}".format(n)), color='red')
			
	return fitPanel
	
def doResPanel(resPanel, rrr, good_rrr, resid, good_resid, maxsma_arcsec, sampling, goodIndexes_integer, badIndexes_integer):
	
	resPanel.set_xlim([maxsma_arcsec*(-0.1),1.1*maxsma_arcsec])
	resPanel.set_ylim([+0.199,-0.199])
	
	minorLocator   = MultipleLocator(xaxisMinorLocator)
	resPanel.xaxis.set_minor_locator(minorLocator)
#  	majorLocator   = MultipleLocator(xaxisMajorLocator)
#  	resPanel.xaxis.set_major_locator(majorLocator)

	majorLocator   = MultipleLocator(0.1)
	minorLocator   = MultipleLocator(0.02)
	resPanel.yaxis.set_major_locator(majorLocator)
	resPanel.yaxis.set_minor_locator(minorLocator)
	resPanel.tick_params(axis = 'y', labelsize=10)
	
 	resPanel.hlines(0, maxsma_arcsec*(-0.1), 1.1*maxsma_arcsec, color='k', linestyles='solid', linewidth=1)
	
	if (sampling == 'log'):
		pylab.plot(rrr, resid, 'wo', markersize=4)
		pylab.plot(good_rrr, good_resid, 'ko', markersize=4)
	elif (sampling == 'comb'):
		good_rrr_undersampled = np.array([])
		good_resid_undersampled = np.array([])
		bad_rrr_undersampled = np.array([])
		bad_resid_undersampled = np.array([])
		for i in range (0,len(rrr),1):
			if i in goodIndexes_integer:
				good_rrr_undersampled = np.append(good_rrr_undersampled, rrr[i])
				good_resid_undersampled = np.append(good_resid_undersampled, resid[i])
			if i in badIndexes_integer:
				bad_rrr_undersampled = np.append(bad_rrr_undersampled, rrr[i])
				bad_resid_undersampled = np.append(bad_resid_undersampled, resid[i])
		resPanel.plot(good_rrr_undersampled, good_resid_undersampled, 'ko', markersize=4)
		resPanel.plot(bad_rrr_undersampled, bad_resid_undersampled, 'wo', markersize=4)
	
	return resPanel
	
def createFigure(bestfitFig, psfFunction, equivalentAxisFit, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, resid, good_resid, componentslist, finalparams, Settings, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing, datatab):
	
	initFitPanel = None
	initResPanel = None
	initEllPanel = None
	#initPAPanel = None
	#initB4Panel = None
	
	if not equivalentAxisFit:
		initFitPanel = plt.subplot2grid((10,11), (0,1), rowspan=3, colspan=5)
		initResPanel = plt.subplot2grid((10,11), (3,1), colspan=5)
		initEllPanel = plt.subplot2grid((10,11), (4,1), rowspan=2, colspan=5)
		plt.setp(initFitPanel.get_xticklabels(), visible=False)
		plt.setp(initResPanel.get_xticklabels(), visible=False)
		initFitPanel.set_ylabel(r'$\mu$ [mag arcsec$^{-2}$]', fontsize=12)	
		initResPanel.set_ylabel(r'$\Delta\mu$', fontsize=14)
		initEllPanel.set_ylabel(r'$\epsilon$', fontsize=18)
		initEllPanel.set_xlabel(r'R$_{\rm maj}$ [arcsec]', fontsize=14)
		
	fitPanel = doFitPanel(initFitPanel, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, componentslist, finalparams, Settings, psfFunction, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing)		
	resPanel = doResPanel(initResPanel, rrr, good_rrr, resid, good_resid, maxsma_arcsec, sampling, goodIndexes_integer, badIndexes_integer)
	ellPanel = doEllipticityPanel(initEllPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	#PAPanel = doPAPanel(initPAPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	#B4Panel = doB4Panel(initB4Panel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	
	#bestfitFig.add_subplot(fitPanel, resPanel, ellPanel, PAPanel, B4Panel)
	bestfitFig.add_subplot(fitPanel, resPanel, ellPanel)
	
	return bestfitFig
	
def doEllipticityPanel(ellPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec):
	
	datatab = datatab[datatab['sma']>2]

	rrr = datatab['sma'] * Settings.pxlToArcsec
	if equivalentAxisFit:
		rrr = datatab['sma'] * Settings.pxlToArcsec * (1 - datatab['ellip'])**0.5

	ellip = datatab['ellip']

	ellPanel.set_xlim([maxsma_arcsec*(-0.1),1.1*maxsma_arcsec])
	minorLocator   = MultipleLocator(xaxisMinorLocator)
	ellPanel.xaxis.set_minor_locator(minorLocator)
#  	majorLocator   = MultipleLocator(xaxisMajorLocator)
#  	ellPanel.xaxis.set_major_locator(majorLocator)
	
	ellPanel.plot(rrr, ellip, color='black', linewidth=2.5)	

	start, end = ellPanel.get_ylim()
	ticks = np.arange(start, end, ((end-start)/6))
	ticks = ticks[1:]
	ellPanel.yaxis.set_ticks(ticks)
	ellPanel.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
	
	return ellPanel
		
def doPAPanel(PAPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec):

	datatab = datatab[datatab['sma']>2]

	rrr = datatab['sma'] * Settings.pxlToArcsec
	if equivalentAxisFit:
		rrr = datatab['sma'] * Settings.pxlToArcsec * (1 - datatab['ellip'])**0.5

	PA = datatab['PA']

	PAPanel.set_xlim([maxsma_arcsec*(-0.1),1.1*maxsma_arcsec])
	minorLocator   = MultipleLocator(xaxisMinorLocator)
	PAPanel.xaxis.set_minor_locator(minorLocator)
#  	majorLocator   = MultipleLocator(xaxisMajorLocator)
#  	PAPanel.xaxis.set_major_locator(majorLocator)
	
	PAPanel.plot(rrr, PA, color='black', linewidth=2.5)	
	
	start, end = PAPanel.get_ylim()
	ticks = np.arange(start, end, ((end-start)/6))
	ticks = ticks[1:]
	PAPanel.yaxis.set_ticks(ticks)
	PAPanel.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

	return PAPanel
		
	
def doB4Panel(B4Panel, datatab, Settings, equivalentAxisFit, maxsma_arcsec):

	datatab = datatab[datatab['sma']>2]

	rrr = datatab['sma'] * Settings.pxlToArcsec
	if equivalentAxisFit:
		rrr = datatab['sma'] * Settings.pxlToArcsec * (1 - datatab['ellip'])**0.5

	B4 = datatab['B4']

	B4Panel.set_xlim([maxsma_arcsec*(-0.1),1.1*maxsma_arcsec])
	minorLocator   = MultipleLocator(xaxisMinorLocator)
	B4Panel.xaxis.set_minor_locator(minorLocator)
# 	majorLocator   = MultipleLocator(xaxisMajorLocator)
# 	B4Panel.xaxis.set_major_locator(majorLocator)

	B4Panel.plot(rrr, B4, color='black', linewidth=2.5)	

	start, end = B4Panel.get_ylim()
	ticks = np.arange(start, end, ((end-start)/6))
	ticks = ticks[1:]
	B4Panel.yaxis.set_ticks(ticks)
	B4Panel.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
	
        xticks = B4Panel.xaxis.get_major_ticks()   #
        xticks[0].label1.set_visible(False)	   #
        xticks[-1].label1.set_visible(False)	   #
        xticks[-2].label1.set_visible(False)	   #
     
        #B4Panel.xaxis.set_ticks(np.arange(0, 800, 200) )

	return B4Panel
	
xaxisMinorLocator = 10
#xaxisMajorLocator = 2000

########## MAIN BODY HERE ############

samplingList = ['comb']

psfList, gaussianSmoothing = createPsf(Settings)

#axisFitList = ['maj', 'eq']
axisFitList = ['maj']

prefixEllipseOutput, galaxy, skyRMS = readDataFileNames('init_n1332.1dfit')
excludedDataFileName = galaxy + '.excl'
excludedRangeList = readExcludedData(excludedDataFileName, Settings) # this is in pixels

skyRMS = float(skyRMS)
skyRMS = 3*skyRMS

bestfitFig = plt.figure()

for sampling in samplingList:	

    #suffix = suffixDict[sampling]										       
    ellipseOutput = prefixEllipseOutput + '_combscale.ell'
    
    for psfFunction in psfList:
    
        for axisFit in axisFitList:	
            
            inputModel = galaxy + '_' + axisFit + '.model'
            
            equivalentAxisFit = False
            if (axisFit == 'eq'):
                equivalentAxisFit = True	
                
            readFitAndPlot(excludedRangeList, galaxy, axisFit, psfFunction, sampling, bestfitFig, equivalentAxisFit)					
    


bestfitFig.subplots_adjust(wspace=0, hspace=0)
bestfitFig.savefig('/Users/gsavorgnan/galaxy_vivisection/papers/ellicular/images/' + galaxy + '_decomposition.eps', format='eps', dpi=1000)
#bestfitFig.savefig(galaxy + '_decomposition.eps', format='eps', dpi=1000)

