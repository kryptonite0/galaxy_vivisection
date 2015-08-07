#### default settings

class Settings:
    
    #performFit = False
    performFit = True
    
    #axisFit = 'maj'
    #axisFit = 'eq'
    #axisFit = 'min'
            
    #fitConvolvedModel = False
    #plotConvolvedFinalModel = False
    #plotConvolvedComponent = False
    fitConvolvedModel = True
    plotConvolvedFinalModel = True
    plotConvolvedComponent = True
    
    useErrorsInFit = False
    #useErrorsInFit = True
    
    smoothing = False
    #smoothing = True
    
    #observation = 'spitzer3.6um'
    #observation = 'n1271' 
    #observation = 'n4342' 
    observation = 'LEDA' 
        
    if (observation == 'n1277'):
        pxlToArcsec = 0.05
        zeropoint = 24.851 # hst for n1277
    
    if (observation == 'n1271'):
        pxlToArcsec = 0.12825
        zeropoint = 24.6949 # Vega
    
    if (observation == 'n4342'):
        pxlToArcsec = 0.046
        zeropoint = 21.639 # Vega f814w
    
    if (observation == 'spitzer3.6um'):
        pxlToArcsec = 1.2232836
        zeropoint = 17.25585 #spitzer 3.6 um

    if (observation == 'LEDA'):
        pxlToArcsec = 0.396
        zeropoint = 22.5 # sdss r

    colordict = {'sersic' : 'red', 'tsersic' : 'red', 'disc' : 'blue', 'gaussian' : 'green', 'psf' : 'purple', 'ferrer' : 'lightblue', 'tdisc' : 'blue', 'gring' : 'grey', 'psfwing' : 'purple'}
    
    #oversamplingStep = 0.2 # pixels
    
#     fitResultsFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/ProFitErole1.0.1/results/1DfitResults.dat'
