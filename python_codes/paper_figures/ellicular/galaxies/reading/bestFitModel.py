from lmfit import Parameters
import numpy as np

from cls.cls import Component

def readBestFitModel(bestFitModelFileName):
    
    componentslist = []
    params = Parameters()
    
    try:
        bestFitModelFile = open(bestFitModelFileName)
        lines = bestFitModelFile.readlines()
        for line in lines:
            if (line.split()[0] != 'Delta' and line.split()[0] != 'Total_galaxy_magnitude'):
                
                comp = Component()
                comp.number = int(line.split()[0])
                comp.name = str(line.split()[1])
                
                #components with 4 parameters
                if (comp.name == 'ferrer'):
                
                    par1name = 'par1_' + str(line.split()[0])
                    par2name = 'par2_' + str(line.split()[0])
                    par3name = 'par3_' + str(line.split()[0])
                    par4name = 'par4_' + str(line.split()[0])
                
                    p1 = (par1name, float(line.split()[2]), None, None, None, None)     # r_out 
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_0 
                    p3 = (par3name, float(line.split()[4]), None, None, None, None)     # alpha
                    p4 = (par4name, float(line.split()[5]), None, None, None, None)     # beta
                            
                    comp.parameters.add_many(p1, p2, p3, p4)  
                    params.add_many(p1, p2, p3, p4)
        
                    componentslist.append(comp)
                
                if (comp.name == 'tsersic'):
                
                    par1name = 'par1_' + str(line.split()[0])
                    par2name = 'par2_' + str(line.split()[0])
                    par3name = 'par3_' + str(line.split()[0])
                    par4name = 'par4_' + str(line.split()[0])
                
                    p1 = (par1name, float(line.split()[2]), None, None, None, None)     # r_e 
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_e 
                    p3 = (par3name, float(line.split()[4]), None, None, None, None)     # n
                    p4 = (par4name, float(line.split()[5]), None, None, None, None)     # r_out
                            
                    comp.parameters.add_many(p1, p2, p3, p4)  
                    params.add_many(p1, p2, p3, p4)
        
                    componentslist.append(comp)
                
                #components with 3 parameters
                if (comp.name == 'sersic'):
                
                    par1name = 'par1_' + str(line.split()[0])
                    par2name = 'par2_' + str(line.split()[0])
                    par3name = 'par3_' + str(line.split()[0])
                
                    p1 = (par1name, float(line.split()[2]), None, None, None, None)     # r_e 
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_e 
                    p3 = (par3name, float(line.split()[4]), None, None, None, None)     # n

                    comp.parameters.add_many(p1, p2, p3)  
                    params.add_many(p1, p2, p3)
        
                    componentslist.append(comp)
                
                if (comp.name == 'tdisc' or comp.name == 'gring'):
                
                    par1name = 'par1_' + str(line.split()[0])
                    par2name = 'par2_' + str(line.split()[0])
                    par3name = 'par3_' + str(line.split()[0])
                
                    p1 = (par1name, float(line.split()[2]), None, None, None, None)     # h # fwhm
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_0 # mu_0 
                    p3 = (par3name, float(line.split()[4]), None, None, None, None)     # r_out # r_0

                    comp.parameters.add_many(p1, p2, p3)  
                    params.add_many(p1, p2, p3)
        
                    componentslist.append(comp)
                
                #components with two parameters	
                elif (comp.name == 'gaussian' or comp.name == 'disc'):
                
                    par1name = 'par1_' + str(line.split()[0])
                    par2name = 'par2_' + str(line.split()[0])
                
                    p1 = (par1name, float(line.split()[2]), None, None, None, None)     # h # fwhm
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_0 # mu_0 

                    comp.parameters.add_many(p1, p2)  
                    params.add_many(p1, p2)
        
                    componentslist.append(comp)
                
                #components with one parameter	
                #elif (comp.name == 'psf' or comp.name == 'psfwing'):
                elif (comp.name == 'psf'):
                
                    par2name = 'par2_' + str(line.split()[0])
                
                    p2 = (par2name, float(line.split()[3]), None, None, None, None)     # mu_0 # mu_0 

                    comp.parameters.add_many(p2)  
                    params.add_many(p2)
        
                    componentslist.append(comp)
               
        bestFitModelFile.close()
        
    except:
        None
        
    return componentslist, params
