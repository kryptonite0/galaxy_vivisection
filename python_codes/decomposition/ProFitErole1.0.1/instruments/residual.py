from instruments.model import buildModel

############# RESIDUAL ##############

def residual(params, x, good_x, y_meas, good_y_meas, good_y_err, goodIndexes, componentslist, psfFunction, fitConvolvedModel, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing):
	
	y_model, good_y_model = buildModel(params, componentslist, x, goodIndexes, psfFunction, fitConvolvedModel, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings, sampling, gaussianSmoothing)
	
	if Settings.useErrorsInFit:
		#res = abs(good_y_meas-good_y_model)/good_y_err
		res = (good_y_meas-good_y_model)**2/good_y_err**2
	else:
		#res = abs(good_y_meas-good_y_model)
		res = (good_y_meas-good_y_model)**2

	return res	
	
	
