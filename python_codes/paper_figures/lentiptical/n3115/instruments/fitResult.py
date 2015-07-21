import numpy as np

from modelComponents.sersic import computeSersicParameters
from modelComponents.disc import computeDiscParameters
from modelComponents.gaussian import computeGaussianParameters
from modelComponents.moffat import computeMoffatParameters

############# PRINT FIT RESULT ##############

def printFitResult(fit, componentslist, psfFunction):
	
	print '--------------------------------------------'
	print 'Results from fit:'
	print '--------------------------------------------'
	print 'Number of iterations:', fit.nfev
	print 'Chisq: ', fit.chisqr
	print 'Reduced chisq: ', fit.redchi
	print 'N_free: ', fit.nfree
	##print 'Delta: ', deltarms   # [sum (residuals**2)/(N-N_dof)]**0.5

	if (fit.success):
		for component in componentslist:
			namepar1 = 'par1_' + str(component.number)
			namepar2 = 'par2_' + str(component.number)
			namepar3 = 'par3_' + str(component.number)
			namepar4 = 'par4_' + str(component.number)
			print 'Component ' + str(component.number) + ': ' + component.name + '.'
			
			if (component.name == 'sersic'):
				r_e = fit.params[namepar1].value
				mu_e = fit.params[namepar2].value
				n = fit.params[namepar3].value				
				b, m_tot = computeSersicParameters(mu_e, r_e, n) 
				print 'R_e [arcsec] = ', ("{0:.2f}".format(r_e)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '	[initial value =', component.parameters[namepar1].value, ']'
				print 'mu_e [mag/arcsec^2] = ', ("{0:.2f}".format(mu_e)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	 [initial value =', component.parameters[namepar2].value, ']'
				print 'n = ', ("{0:.2f}".format(n)), '+/-', ("{0:.2f}".format(fit.params[namepar3].stderr)), '       [initial value =', component.parameters[namepar3].value, ']'
				print 'b_n = ', ("{0:.2f}".format(b))
				print 'm_tot [mag] = ', ("{0:.2f}".format(m_tot))
			if (component.name == 'disc'):
				h = fit.params[namepar1].value
				mu_0 = fit.params[namepar2].value
				m_tot = computeDiscParameters(mu_0, h) 
				print 'h [arcsec] = ', ("{0:.2f}".format(h)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '       [initial value =', component.parameters[namepar1].value, ']'
	 			print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	 [initial value =', component.parameters[namepar2].value, ']'
				print 'm_tot [mag] = ', ("{0:.2f}".format(m_tot))
			if (component.name == 'gaussian'):
				fwhm = fit.params[namepar1].value
				mu_0 = fit.params[namepar2].value
				m_tot = computeGaussianParameters(fwhm, mu_0) 
				print 'fwhm [arcsec] = ', ("{0:.2f}".format(fwhm)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '	 [initial value =', component.parameters[namepar1].value, ']'
	 			print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	[initial value =', component.parameters[namepar2].value, ']'
				print 'm_tot [mag] = ', ("{0:.2f}".format(m_tot))
			if (component.name == 'psf'):
	 			mu_0 = fit.params[namepar2].value
				if (psfFunction.name == 'gaussian'):
					m_tot = computeGaussianParameters(psfFunction.gaussianFWHM, mu_0) 
				elif (psfFunction.name == 'moffat'):
					alpha = psfFunction.moffatAlpha
					beta = psfFunction.moffatBeta
					m_tot = computeMoffatParameters(alpha, beta, mu_0) 
				print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	[initial value =', component.parameters[namepar2].value, ']'
	 			print 'm_tot [mag] = ', ("{0:.2f}".format(m_tot))
			if (component.name == 'ferrer'):
				r_out = fit.params[namepar1].value
				mu_0 = fit.params[namepar2].value
				alpha = fit.params[namepar3].value
				beta = fit.params[namepar4].value
				print 'R_out [arcsec] = ', ("{0:.2f}".format(r_out)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '	[initial value =', component.parameters[namepar1].value, ']'
				print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	 [initial value =', component.parameters[namepar2].value, ']'
				print 'alpha [pixel] = ', ("{0:.2f}".format(alpha)), '+/-', ("{0:.2f}".format(fit.params[namepar3].stderr)), '       [initial value =', component.parameters[namepar3].value, ']'
				print 'beta = ', ("{0:.2f}".format(beta)), '+/-', ("{0:.2f}".format(fit.params[namepar4].stderr)), '       [initial value =', component.parameters[namepar4].value, ']'
				print 'm_tot [mag] = -99'
			if (component.name == 'tdisc'):
				h = fit.params[namepar1].value
				mu_0 = fit.params[namepar2].value
				r_out = fit.params[namepar3].value
				print 'h [arcsec] = ', ("{0:.2f}".format(h)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '       [initial value =', component.parameters[namepar1].value, ']'
	 			print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', fit.params[namepar2].stderr, '	 [initial value =', component.parameters[namepar2].value, ']'
				print 'R_out [arcsec] = ', ("{0:.2f}".format(r_out)), '+/-', ("{0:.2f}".format(fit.params[namepar3].stderr)), '	 [initial value =', component.parameters[namepar3].value, ']'
				print 'm_tot [mag] = 99'
			if (component.name == 'tsersic'):
				r_e = fit.params[namepar1].value
				mu_e = fit.params[namepar2].value
				n = fit.params[namepar3].value
				r_out = fit.params[namepar4].value
				b, m_tot = computeSersicParameters(mu_e, r_e, n) 
				print 'R_e [arcsec] = ', ("{0:.2f}".format(r_e)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '	[initial value =', component.parameters[namepar1].value, ']'
				print 'mu_e [mag/arcsec^2] = ', ("{0:.2f}".format(mu_e)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	 [initial value =', component.parameters[namepar2].value, ']'
				print 'n = ', ("{0:.2f}".format(n)), '+/-', ("{0:.2f}".format(fit.params[namepar3].stderr)), '       [initial value =', component.parameters[namepar3].value, ']'
				print 'R_out [arcsec] = ', ("{0:.2f}".format(r_out)), '+/-', ("{0:.2f}".format(fit.params[namepar4].stderr)), '	[initial value =', component.parameters[namepar4].value, ']'
				print 'b_n = ', ("{0:.2f}".format(b))
				print 'm_tot [mag] = 99'
			if (component.name == 'gring'):
				fwhm = fit.params[namepar1].value
				mu_0 = fit.params[namepar2].value
				r_0 = fit.params[namepar3].value
				print 'fwhm [arcsec] = ', ("{0:.2f}".format(fwhm)), '+/-', ("{0:.2f}".format(fit.params[namepar1].stderr)), '	 [initial value =', component.parameters[namepar1].value, ']'
	 			print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	[initial value =', component.parameters[namepar2].value, ']'
				print 'R_0 [arcsec] = ', ("{0:.2f}".format(r_0)), '+/-', ("{0:.2f}".format(fit.params[namepar3].stderr)), '	[initial value =', component.parameters[namepar3].value, ']'
				print 'm_tot [mag] = -99'
			if (component.name == 'psfwing'):
	 			mu_0 = fit.params[namepar2].value
				print 'mu_0 [mag/arcsec^2] = ', ("{0:.2f}".format(mu_0)), '+/-', ("{0:.2f}".format(fit.params[namepar2].stderr)), '	[initial value =', component.parameters[namepar2].value, ']'
	 			print 'm_tot [mag] = -99'
			print	
	
		print '--------------------------------------------'
		print '--------------------------------------------'
		print
	else:
		print 'Unsuccesful fit :('
		print '--------------------------------------------'
		print '--------------------------------------------'
		print
		
