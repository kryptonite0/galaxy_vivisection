import numpy as np

######## READ ELLIPSE OUTPUT IN ARRAYS #############

def readEllipseOutput(data):

      data = open(data)
      
      datatab = np.genfromtxt(data, 
              dtype={'names': ('row', 'sma', 'intens', 'intens_err', 'pix_var', 'rms', 'ellip', 'ellip_err', 'PA', 'PA_err', 'X0', 'X0_ERR', 'Y0', 'Y0_ERR', 'GRAD', 'GRAD_ERR', 'GRAD_R_ERR', 'RSMA', 'MAG', 'MAG_LERR', 'MAG_UERR', 'TFLUX_E', 'TFLUX_C', 'TMAG_E', 'TMAG_C', 'NPIX_E', 'NPIX_C', 'A3', 'A3_ERR', 'B3', 'B3_ERR', 'A4', 'A4_ERR', 'B4', 'B4_ERR'),
              'formats': ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, skiprows=5, missing_values=('INDEF'), filling_values=0.0)      
      # !!! add missing value control (see genfromtxt function)
      data.close()
      
      return datatab

