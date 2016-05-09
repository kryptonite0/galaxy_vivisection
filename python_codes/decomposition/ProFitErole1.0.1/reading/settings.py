def readSettings(txt):

    data = open(txt)
    
    for line in data:
        prefixEllipseOutput = str(line.split()[0])
        inputModel = str(line.split()[1])
        skyRMS = float(line.split()[2])	
        
    return prefixEllipseOutput, inputModel, skyRMS	
