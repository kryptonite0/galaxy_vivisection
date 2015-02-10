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


import numpy as np
import pylab
from matplotlib.pyplot import fill_between
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc, font_manager

from textwrap import wrap


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
		fitPanel.text(maxsma_arcsec*0.7, (minmu+0.5), r'$R_{\rm e} = $' + ("{0:.1f}".format(r_e)), color='red')
		fitPanel.text(maxsma_arcsec*0.7, (minmu+1.5), r'$\mu_{\rm e} = $' + ("{0:.2f}".format(mu_e)), color='red')
		fitPanel.text(maxsma_arcsec*0.7, (minmu+2.5), r'$n = $' + ("{0:.1f}".format(n)), color='red')
			
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
	
def makePaperFigure(bestfitFig, psfFunction, equivalentAxisFit, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, resid, good_resid, componentslist, finalparams, Settings, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing, datatab):
	
	initFitPanel = None
	initResPanel = None
	initEllPanel = None
	initPAPanel = None
	initB4Panel = None
	
	if not equivalentAxisFit:
		initFitPanel = plt.subplot2grid((10,11), (0,1), rowspan=3, colspan=5)
		initResPanel = plt.subplot2grid((10,11), (3,1), colspan=5)
		initEllPanel = plt.subplot2grid((10,11), (4,1), rowspan=2, colspan=5)
		initPAPanel = plt.subplot2grid((10,11), (6,1), rowspan=2, colspan=5)
		initB4Panel = plt.subplot2grid((10,11), (8,1), rowspan=2, colspan=5)
		plt.setp(initFitPanel.get_xticklabels(), visible=False)
		plt.setp(initResPanel.get_xticklabels(), visible=False)
		plt.setp(initEllPanel.get_xticklabels(), visible=False)
		plt.setp(initPAPanel.get_xticklabels(), visible=False)
		initFitPanel.set_ylabel(r'$\mu$ [mag arcsec$^{-2}$]')	
		initResPanel.set_ylabel(r'$\Delta\mu$')
		initEllPanel.set_ylabel(r'$\epsilon$')
		initPAPanel.set_ylabel(r'$PA$ [deg]')
		initB4Panel.set_ylabel(r'$B4$')
		initB4Panel.set_xlabel(r'R$_{\rm maj}$ [arcsec]', fontsize=14)
	if equivalentAxisFit:	
		initFitPanel = plt.subplot2grid((10,11), (0,6), rowspan=3, colspan=5)
		initResPanel = plt.subplot2grid((10,11), (3,6), colspan=5)
		initEllPanel = plt.subplot2grid((10,11), (4,6), rowspan=2, colspan=5)
		initPAPanel = plt.subplot2grid((10,11), (6,6), rowspan=2, colspan=5)
		initB4Panel = plt.subplot2grid((10,11), (8,6), rowspan=2, colspan=5)
		plt.setp(initFitPanel.get_xticklabels(), visible=False)
		plt.setp(initResPanel.get_xticklabels(), visible=False)
		plt.setp(initEllPanel.get_xticklabels(), visible=False)
		plt.setp(initPAPanel.get_xticklabels(), visible=False)		
		plt.setp(initFitPanel.get_yticklabels(), visible=False)
		plt.setp(initResPanel.get_yticklabels(), visible=False)
		plt.setp(initEllPanel.get_yticklabels(), visible=False)
		plt.setp(initPAPanel.get_yticklabels(), visible=False)
		plt.setp(initB4Panel.get_yticklabels(), visible=False)
		initB4Panel.set_xlabel(r'R$_{\rm eq}$ [arcsec]', fontsize=14)
			
	fitPanel = doFitPanel(initFitPanel, rrr, mu, good_rrr, good_mu, bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, componentslist, finalparams, Settings, psfFunction, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing)		
	resPanel = doResPanel(initResPanel, rrr, good_rrr, resid, good_resid, maxsma_arcsec, sampling, goodIndexes_integer, badIndexes_integer)
	ellPanel = doEllipticityPanel(initEllPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	PAPanel = doPAPanel(initPAPanel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	B4Panel = doB4Panel(initB4Panel, datatab, Settings, equivalentAxisFit, maxsma_arcsec) 
	
	bestfitFig.add_subplot(fitPanel, resPanel, ellPanel, PAPanel, B4Panel)
	
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
	
	return B4Panel
	
xaxisMinorLocator = 200
xaxisMajorLocator = 2000