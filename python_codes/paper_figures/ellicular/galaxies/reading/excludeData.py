def readExcludedData(txt, Settings):
    
    excludedRangeList = None
    
    try:
        
        excludedRangeList = [] 
        data = open(txt)
    
        for line in data:
            
            pt1 = float(line.split()[0]) * Settings.pxlToArcsec
            pt2 = float(line.split()[1]) * Settings.pxlToArcsec
            excludedRange = [pt1, pt2]	
            excludedRangeList.append(excludedRange)
    except:
        None
        
    return excludedRangeList
    
