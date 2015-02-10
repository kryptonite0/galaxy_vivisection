import numpy as np

############# BUILD PSF WING MODEL ##############

def buildPsfWing(x, params, component, sampling, psfwing_02pxscale_datatab, psfwing_logscale_datatab, Settings):
  
  if (sampling == 'comb'):	     		
    namepar2 = 'par2_' + str(component.number) # mu_0
#     x_e = x[x>max(psfwing_02pxscale_datatab['sma'])]
#     y_psfwing_i = psfwing_02pxscale_datatab['intens'][(psfwing_02pxscale_datatab['sma'])>=min(x)]
    x_e = x[len(psfwing_02pxscale_datatab['sma']):]
    y_psfwing_i = psfwing_02pxscale_datatab['intens'][:len(psfwing_02pxscale_datatab['sma'])]

  elif (sampling == 'log'):
    namepar2 = 'par2_' + str(component.number) # mu_0 
#     x_e = x[x>max(psfwing_logscale_datatab['sma'])]
#     y_psfwing_i = psfwing_logscale_datatab['intens'][(psfwing_logscale_datatab['sma'])>=min(x)]
    x_e = x[len(psfwing_logscale_datatab['sma']):]
    y_psfwing_i = psfwing_logscale_datatab['intens'][:len(psfwing_logscale_datatab['sma'])]
  
  y_psfwing_e = [0.0] * x_e		
  y_psfwing = np.concatenate((y_psfwing_i, y_psfwing_e))
  y_psfwing = 10**(Settings.zeropoint/2.5) * (10**(-params[namepar2].value/2.5)) * y_psfwing
  
  #print 'quiiiii', len(y_psfwing)
  
  return y_psfwing


