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
    observation = 'm1216' 
    #observation = 'n1271'
        
    if (observation == 'n1277'):
        pxlToArcsec = 0.05
        zeropoint = 24.851 # hst for n1277
    
    if (observation == 'n1271'):
        pxlToArcsec = 0.12825
        zeropoint = 24.6949 # Vega
    
    if (observation == 'm1216'):
        pxlToArcsec = 0.12825
        zeropoint = 24.6949 # Vega
    
    if (observation == 'spitzer3.6um'):
        pxlToArcsec = 1.2232836
        zeropoint = 17.25585 #spitzer 3.6 um

    colordict = {'sersic' : 'red', 'tsersic' : 'red', 'disc' : 'blue', 'gaussian' : 'green', 'psf' : 'purple', 'ferrer' : 'lightblue', 'tdisc' : 'blue', 'gring' : 'green', 'psfwing' : 'purple', 'sersicdisc' : 'blue'}
    linestyledict = {'sersic' : '-', 'tsersic' : '-', 'disc' : '--', 'gaussian' : '-', 'psf' : '-', 'ferrer' : '-', 'tdisc' : '-', 'gring' : ':', 'psfwing' : '-', 'sersicdisc' : '--'}
    
    #oversamplingStep = 0.2 # pixels
    
#     fitResultsFileName = '/Users/gsavorgnan/Dropbox/giulia_e_basta/ProFitErole1.0.1/results/1DfitResults.dat'
