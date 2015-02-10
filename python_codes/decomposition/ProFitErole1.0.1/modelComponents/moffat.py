import numpy as np

def computeMoffatParameters(alpha, beta, mu_0):
    
    m_tot = mu_0 - 2.5*np.log10(np.pi*(alpha)**2/(beta-1))
    
    return m_tot