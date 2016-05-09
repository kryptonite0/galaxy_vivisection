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

from reading.dataFileNames import readInitialSettings
from reading.ellipseOutput import readEllipseOutput
from reading.inputModel import readInputModel
from reading.excludeData import readExcludedData

from cls.cls import *
from cls.settings import *
from cls.psf import *

from instruments.fitResult import printFitResult
from instruments.model import buildModel
from instruments.output import *

from modelComponents.disc import buildDisc
from modelComponents.ferrer import buildFerrer
from modelComponents.gaussian import buildGaussian
from modelComponents.psf import buildPsf
from modelComponents.psfWing import buildPsfWing
from modelComponents.sersic import buildSersic
from modelComponents.truncDisc import buildTruncDisc
from modelComponents.truncSersic import buildTruncSersic
from modelComponents.gaussRing import buildGaussRing

from instruments.residual import residual
from instruments import b_n
from instruments.plotting import *

######### LIMITS OF FIT PARAMETERS #############

par1min = 0.0000001   # r_e or h or fwhm ...
par1max = None   

par2min = None        # mu_e or mu_0 or ...
par2max = 35.0   

par3min = 0.001       # n, alphaferrer
par3max = 20.0

par4min = 0.001       # betaferrer
par4max = 4.0


####################################################

terminal = sys.stdout


########## FIT AND PLOT ################

def performFitAndPlot(excludedRangeList, galaxyName, axisFit, psfFunction, sampling, bestfitFig, equivalentAxisFit):

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

	componentslist, params, psfwing_02pxscale_datatab, psfwing_logscale_datatab = readInputModel(inputModel, 
		equivalentAxisFit, Settings)
	
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
# 	if (minorAxisFit):
# 		good_rrr = good_rrr * (1. - good_ellip) # semi minor axis	
# 		bad_rrr = bad_rrr * (1. - bad_ellip)
# 		rrr = rrr *(1. - ellip) 
		
	# perform minimization

	if (Settings.performFit):
		
		# perform minimization (lm algorithm used)
		fit = minimize(residual, params, args=(rrr, good_rrr, mu, good_mu, good_mu_up_err, goodIndexes, 
			componentslist, psfFunction, Settings.fitConvolvedModel, psfwing_02pxscale_datatab, 
			psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing), method='leastsq')
		# report fit
		#printFitResult(fit, componentslist, psfFunction)
		# calculate final result
		finalparams = fit.params
		
	else:
		finalparams = params
	
	
	y_finalmodel_sb, good_y_finalmodel_sb = buildModel(finalparams, componentslist, rrr, goodIndexes, 
		psfFunction, Settings.plotConvolvedFinalModel, psfwing_02pxscale_datatab, 
		psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing)
	resid = mu - y_finalmodel_sb
	good_resid = good_mu - good_y_finalmodel_sb

	#compute delta rms
	deltarms = 0
	for i in range(0,len(good_rrr),1):
		deltarms = deltarms + (good_resid[i])**2
	deltarms = (deltarms / fit.nfree)**0.5	
	
	
			
	#if happyOfFit and Settings.performFit:		
	if Settings.performFit:		
		logfile = open('1dfit.log', 'a+')
		sys.stdout = logfile
		print 'Input ellipse file: ', ellipseOutput
		print 'Excluded data: ', excludedRangeList
		print 'Axis Used: ', axisFit
		print 'x-Axis units: arcsec'
		print 'pixels to arcsec scale: ', Settings.pxlToArcsec
		print 'Smoothing: ', Settings.smoothing, 'sigma = ', (gaussianSmoothing.gaussianFWHM / (2.3548 * Settings.pxlToArcsec))
		print 'Zeropoint: ', Settings.zeropoint 
		print 'Sky RMS: ', skyRMS
		print 'Delta: ', deltarms   # [sum (residuals**2)/(N-N_dof)]**0.5
		print 'Convolution done? ', Settings.fitConvolvedModel 
		print 'PSF: ', psfFunction.name
		print 'Fit weighted with errors? ', Settings.useErrorsInFit 	
					
		printFitResult(fit, componentslist, psfFunction)
		
		logfile.close()
		
		sys.stdout = terminal
		print 'Input ellipse file: ', ellipseOutput
		print 'Excluded data: ', excludedRangeList
		print 'Axis Used: ', axisFit
		print 'x-Axis units: arcsec'
		print 'pixels to arcsec scale: ', Settings.pxlToArcsec
		print 'Smoothing: ', Settings.smoothing, 'sigma = ', (gaussianSmoothing.gaussianFWHM / (2.3548 * Settings.pxlToArcsec))
		print 'Zeropoint: ', Settings.zeropoint 
		print 'Sky RMS: ', skyRMS
		print 'Delta: ', deltarms   # [sum (residuals**2)/(N-N_dof)]**0.5
		print 'Convolution done? ', Settings.fitConvolvedModel 
		print 'PSF: ', psfFunction.name
		print 'Fit weighted with errors? ', Settings.useErrorsInFit 	
			
		printFitResult(fit, componentslist, psfFunction)
		
		
		produceOutputFile(Settings, galaxyName, axisFit, psfFunction, sampling, componentslist, finalparams, deltarms)
		sys.stdout = terminal
		
#### produce figures here	
	
	#bestfitFig = addFitResPanels(bestfitFig, psfFunction, equivalentAxisFit, rrr, mu, good_rrr, good_mu, 
		#bad_rrr, bad_mu, maxsma_arcsec, minmu, maxmu, resid, good_resid, componentslist, finalparams, 
		#Settings, psfwing_02pxscale_datatab, psfwing_logscale_datatab, goodIndexes, skyRMS, deltarms, 
		#sampling, goodIndexes_integer, badIndexes_integer, gaussianSmoothing)
					
	
########## MAIN BODY HERE ############

#samplingList = ['lin', 'log']
samplingList = ['lin']

suffixDict = {'lin' : '_linscale', 'log' : '_logscale'}

axisFitList = ['maj', 'eq']
#axisFitList = ['maj']
#axisFitList = ['eq']

galaxyName, prefixEllipseOutput, skyRMS, Settings, psfFunction, gaussianSmoothing = readInitialSettings('settings.1dfit', Settings)

excludedDataFileName = galaxyName + '.excl'
excludedRangeList = readExcludedData(excludedDataFileName, Settings) # this is in pixels

skyRMS = float(skyRMS)
skyRMS = 3*skyRMS

bestfitFig = plt.figure()

for sampling in samplingList:	

	suffix = suffixDict[sampling]										       
	ellipseOutput = prefixEllipseOutput + suffix + '.ell'
	
	#for psfFunction in psfList:
	
	for axisFit in axisFitList:	
		
		inputModel = galaxyName + '_' + axisFit + '.model'
		
		equivalentAxisFit = False
		if (axisFit == 'eq'):
			equivalentAxisFit = True	
			
		performFitAndPlot(excludedRangeList, galaxyName, axisFit, psfFunction, sampling, bestfitFig, 
			equivalentAxisFit)					
	
	#bestfitFig.subplots_adjust(wspace=0, hspace=0)
	#bestfitFig.savefig(galaxy + 'bestfit_' + sampling + '.eps', format='eps', dpi=1000)
	
#ellipseOutput = prefixEllipseOutput + '_combscale.ell'
#datatab = readEllipseOutput(ellipseOutput) 

#exclude R=0 point
#datatab = datatab[datatab['sma']>0]
#exclude data fainter than (3) sky RMS
#datatab = datatab[datatab['intens']>skyRMS]

#isophAnalFig = plt.figure()

#isophAnalFig = addPanelsToIsophotalAnalysisFigure(isophAnalFig, datatab, Settings)
#isophAnalFig.savefig(galaxy + 'isophotalAnalysis.eps', format='eps', dpi=1000)

