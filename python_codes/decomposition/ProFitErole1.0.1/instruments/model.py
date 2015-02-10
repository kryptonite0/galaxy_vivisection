import numpy as np

from modelComponents.disc import buildDisc
from modelComponents.ferrer import buildFerrer
from modelComponents.gaussian import buildGaussian
from modelComponents.psf import buildPsf
from modelComponents.psfWing import buildPsfWing
from modelComponents.sersic import buildSersic
from modelComponents.truncDisc import buildTruncDisc
from modelComponents.truncSersic import buildTruncSersic
from modelComponents.gaussRing import buildGaussRing

from convolution.psfConvolution import doPsfConvolution

############# BUILD MODEL ##############

def buildModel(params, componentslist, xxx, goodIndexes, psfFunction, convolve, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing):
						
	model = 0
	for component in componentslist:
		if (component.name == 'sersic'):
			y_component = buildSersic(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'disc'):
			y_component = buildDisc(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'gaussian'):
			y_component = buildGaussian(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'psf' or component.name == 'psfwing'):	
			y_component = xxx*0
		elif (component.name == 'ferrer'):	
			y_component = buildFerrer(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'tdisc'):	
			y_component = buildTruncDisc(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'tsersic'):	
			y_component = buildTruncSersic(xxx, params, component, psfFunction, False, Settings.zeropoint)
		elif (component.name == 'gring'):
			y_component = buildGaussRing(xxx, params, component, psfFunction, False, Settings.zeropoint)		
		else:
			y_component = None
			print 'ERROR: Not valid function name: ', component.name
		model = model + y_component
	
	if convolve:
		model = doPsfConvolution(model, psfFunction, xxx)		
		if Settings.smoothing:
			model = doPsfConvolution(model, gaussianSmoothing, xxx)
	
	
	for component in componentslist:
		if (component.name == 'psf'):
			y_psf = buildPsf(xxx, params, component, psfFunction, Settings.zeropoint)
			model = model + y_psf
		if (component.name == 'psfwing'):
			y_psfwing = buildPsfWing(xxx, params, component, sampling, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings)
			model = model + y_psfwing	
	
	model_sb = Settings.zeropoint - 2.5*np.log10(model)
		
	good_model_sb = np.array([])
	
	for i in range(0,len(xxx),1):
		if i in goodIndexes:
			good_model_sb = np.append(good_model_sb, model_sb[i])		
	
	return model_sb, good_model_sb

