import numpy as np
import sys
from scipy.special import gamma
import os.path

from modelComponents.sersic import computeSersicParameters
from modelComponents.disc import computeDiscParameters
from modelComponents.gaussian import computeGaussianParameters
from modelComponents.moffat import computeMoffatParameters
from modelComponents.ferrer import computeFerrerParameters
from modelComponents.gaussRing import computeGaussRingParameters

def produceOutputFile(Settings, galaxy, axisFit, psfFunction, sampling, componentslist, finalparams, deltarms):
    
    outputFileName = galaxy + '_' + axisFit + '_' + psfFunction.name + '_' + sampling + '_par_SM.dat'
    outputFile = open(outputFileName, 'w')
    sys.stdout = outputFile
    
    flux_total_galaxy = 0
    
    for component in componentslist:
        namepar1 = 'par1_' + str(component.number)
        namepar2 = 'par2_' + str(component.number)
        namepar3 = 'par3_' + str(component.number)
        namepar4 = 'par4_' + str(component.number)
            
        if (component.name == 'sersic'):
            r_e = finalparams[namepar1].value
            mu_e = finalparams[namepar2].value
            n = finalparams[namepar3].value
            b, m_tot = computeSersicParameters(mu_e, r_e, n)
            print component.number,
            print component.name,
            print("{0:.4f}".format(r_e)),
            print("{0:.4f}".format(mu_e)),
            print("{0:.4f}".format(n)),
            print '-99.9999',
            print("{0:.4f}".format(m_tot))
        if (component.name == 'disc'):
            h= finalparams[namepar1].value
            mu_0 = finalparams[namepar2].value
            m_tot = computeDiscParameters(mu_0, h) 
            print component.number,
            print component.name,
            print("{0:.4f}".format(h)),
            print("{0:.4f}".format(mu_0)),
            print '-99.9999',
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'gaussian'):
            fwhm = finalparams[namepar1].value
            mu_0 = finalparams[namepar2].value
            m_tot = computeGaussianParameters(fwhm, mu_0)
            print component.number,
            print component.name,
            print("{0:.4f}".format(fwhm)),
            print("{0:.4f}".format(mu_0)),
            print '-99.9999',
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'psf'):
            mu_0 = finalparams[namepar2].value
            if (psfFunction.name == 'gaussian'):
                m_tot = computeGaussianParameters(psfFunction.gaussianFWHM, mu_0)
            elif (psfFunction.name == 'moffat'):
                alpha = psfFunction.moffatAlpha
                beta = psfFunction.moffatBeta
                m_tot = computeMoffatParameters(alpha, beta, mu_0)
            print component.number,
            print component.name,
            print '-99.9999',          
            print("{0:.4f}".format(mu_0)),
            print '-99.9999',
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'ferrer'):
            r_out = finalparams[namepar1].value
            mu_0 = finalparams[namepar2].value
            alpha = finalparams[namepar3].value
            beta = finalparams[namepar4].value
            m_tot = computeFerrerParameters(r_out, mu_0, alpha, beta)
            print component.number,
            print component.name,
            print("{0:.4f}".format(r_out)),
            print("{0:.4f}".format(mu_0)),
            print("{0:.4f}".format(alpha)),
            print("{0:.4f}".format(beta)),
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'tdisc'):
            h = finalparams[namepar1].value
            mu_0 = finalparams[namepar2].value
            r_out = finalparams[namepar3].value
            m_tot = -99.9999
            print component.number,
            print component.name,
            print("{0:.4f}".format(h)),
            print("{0:.4f}".format(mu_0)),
            print("{0:.4f}".format(r_out)),
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'tsersic'):
            r_e = finalparams[namepar1].value
            mu_e = finalparams[namepar2].value
            n = finalparams[namepar3].value
            r_out = finalparams[namepar4].value
            m_tot = -99.9999
            print component.number,
            print component.name,
            print("{0:.4f}".format(r_e)),
            print("{0:.4f}".format(mu_e)),
            print("{0:.4f}".format(n)),
            print("{0:.4f}".format(r_out)),
            print("{0:.4f}".format(m_tot))
        if (component.name == 'gring'):
            fwhm = finalparams[namepar1].value
            mu_0 = finalparams[namepar2].value
            r_0 = finalparams[namepar3].value
            m_tot = computeGaussRingParameters(fwhm, mu_0, r_0)
            print component.number,
            print component.name,
            print("{0:.4f}".format(fwhm)),
            print("{0:.4f}".format(mu_0)),
            print("{0:.4f}".format(r_0)),
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        if (component.name == 'psfwing'):
            mu_0 = finalparams[namepar2].value
            print component.number,
            print component.name,
            m_tot = -99.9999
            print '-99.9999',
            print("{0:.2f}".format(mu_0)),
            print '-99.9999',
            print '-99.9999',
            print("{0:.4f}".format(m_tot)) 
        
        flux_total_galaxy = flux_total_galaxy + 10**(-0.4*m_tot)

#     try:
#         if (componentslist[0].number == 1 and componentslist[0].name == 'sersic'):
#             r_e = finalparams['par1_1'].value
#             mu_e = finalparams['par2_1'].value
#             n = finalparams['par3_1'].value
#             b, m_tot = computeSersicParameters(mu_e, r_e, n)
#             print 'sersic',
#             print("{0:.4f}".format(r_e)),
#             print("{0:.4f}".format(mu_e)),
#             print("{0:.4f}".format(n)),
#             print("{0:.4f}".format(b)),
#             print("{0:.4f}".format(m_tot))
#     except:
#         print 'Sersic function has not been fit'        
#     try:
#         if (componentslist[1].number == 2 and componentslist[1].name == 'disc'):
#             h = finalparams['par1_2'].value
#             mu_0 = finalparams['par2_2'].value
#             m_tot = computeDiscParameters(mu_0, h) 
#             print 'disc',
#             print("{0:.4f}".format(h)),
#             print("{0:.4f}".format(mu_0)),
#             print("{0:.4f}".format(m_tot)) 
#     except:
#         print 'Disc function has not been fit'        
    
    print 'Delta', deltarms    
    print 'Total_galaxy_magnitude', -2.5*np.log10(flux_total_galaxy) 
