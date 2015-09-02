from lmfit import Parameters
import numpy as np

from cls.cls import *

from reading.ellipseOutput import readEllipseOutput

def readInputModel(txt, equivalentAxisFit, Settings):

	psfwing_02pxscale_datatab = None
	psfwing_logscale_datatab = None

	componentslist = []
	params = Parameters()
	
	data = open(txt)
	for line in data:
		if (line[0] != '#'):
			comp = Component()
			comp.number = int(line.split()[0])
			comp.name = str(line.split()[1])
			
			#components with 4 parameters
			if (comp.name == 'ferrer'):
			
				par1name = 'par1_' + str(line.split()[0])
				par2name = 'par2_' + str(line.split()[0])
				par3name = 'par3_' + str(line.split()[0])
				par4name = 'par4_' + str(line.split()[0])
			
				p1 = (par1name, float(line.split()[2]), True, 0.01, None, None)     # r_out 
				if (line.split()[3] == 'False'):
					p1 = (par1name, float(line.split()[2]), False, 0.01, None, None) 
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)     # mu_0 
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				p3 = (par3name, float(line.split()[6]), True, 0.01, 4.0, None)     # alpha
				if (line.split()[7] == 'False'):
					p3 = (par3name, float(line.split()[6]), False, 0.01, 4.0, None) 	
				p4 = (par4name, float(line.split()[8]), True, 0.01, 1.999, None)     # beta
				if (line.split()[9] == 'False'):
					p4 = (par4name, float(line.split()[8]), False, 0.01, 1.999, None) 	
						
				comp.parameters.add_many(p1, p2, p3, p4)  
				params.add_many(p1, p2, p3, p4)
	
				componentslist.append(comp)
			
			if (comp.name == 'tsersic'):
			
				par1name = 'par1_' + str(line.split()[0])
				par2name = 'par2_' + str(line.split()[0])
				par3name = 'par3_' + str(line.split()[0])
				par4name = 'par4_' + str(line.split()[0])
			
				p1 = (par1name, float(line.split()[2]), True, 0.01, None, None)           # r_e 
				if (line.split()[3] == 'False'):
					p1 = (par1name, float(line.split()[2]), False, 0.01, None, None) 
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)           # mu_e 
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				p3 = (par3name, float(line.split()[6]), True, 0.01, 20.0, None)           # n
				if (line.split()[7] == 'False'):
					p3 = (par3name, float(line.split()[6]), False, 0.01, 20.0, None) 		
				p4 = (par4name, float(line.split()[8]), True, 0.01, None, None)           # r_out 
				if (line.split()[9] == 'False'):
					p4 = (par4name, float(line.split()[8]), False, 0.01, None, None) 
						
				comp.parameters.add_many(p1, p2, p3, p4)  
				params.add_many(p1, p2, p3, p4)
	
				componentslist.append(comp)
			
			#components with 3 parameters
			if (comp.name == 'sersic' or comp.name == 'sersicdisc'):
			
				par1name = 'par1_' + str(line.split()[0])
				par2name = 'par2_' + str(line.split()[0])
				par3name = 'par3_' + str(line.split()[0])
			
				p1 = (par1name, float(line.split()[2]), True, 0.01, None, None)     # r_e 
				if (line.split()[3] == 'False'):
					p1 = (par1name, float(line.split()[2]), False, 0.01, None, None) 
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)     # mu_e 
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				p3 = (par3name, float(line.split()[6]), True, 0.01, 20.0, None)     # n
				if (line.split()[7] == 'False'):
					p3 = (par3name, float(line.split()[6]), False, 0.01, 20.0, None) 		
				comp.parameters.add_many(p1, p2, p3)  
				params.add_many(p1, p2, p3)
	
				componentslist.append(comp)
			
			if (comp.name == 'tdisc' or comp.name == 'gring'):
			
				par1name = 'par1_' + str(line.split()[0])
				par2name = 'par2_' + str(line.split()[0])
				par3name = 'par3_' + str(line.split()[0])
			
				p1 = (par1name, float(line.split()[2]), True, 0.01, None, None)     # h # fwhm
				if (line.split()[3] == 'False'):
					p1 = (par1name, float(line.split()[2]), False, 0.01, None, None) 
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)     # mu_0 # mu_0 
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				p3 = (par3name, float(line.split()[6]), True, 0.01, None, None)     # r_out # r_0
				if (line.split()[7] == 'False'):
					p3 = (par3name, float(line.split()[6]), False, 0.01, None, None) 		
				comp.parameters.add_many(p1, p2, p3)  
				params.add_many(p1, p2, p3)
	
				componentslist.append(comp)
			
			#components with two parameters	
			elif (comp.name == 'gaussian' or comp.name == 'disc'):
			
				par1name = 'par1_' + str(line.split()[0])
				par2name = 'par2_' + str(line.split()[0])
			
				p1 = (par1name, float(line.split()[2]), True, 0.01, None, None)     # h or fwhm ..
				if (line.split()[3] == 'False'):
					p1 = (par1name, float(line.split()[2]), False, 0.01, None, None) 
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)     # mu_0 or mag
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				comp.parameters.add_many(p1, p2)  
				params.add_many(p1, p2)
	
				componentslist.append(comp)
			
			#components with one parameter	
			elif (comp.name == 'psf' or comp.name == 'psfwing'):
			
				par2name = 'par2_' + str(line.split()[0])
			
				p2 = (par2name, float(line.split()[4]), True, None, 35.0, None)     # mu_0 or mag
				if (line.split()[5] == 'False'):
					p2 = (par2name, float(line.split()[4]), False, None, 35.0, None)    		
				comp.parameters.add_many(p2)  
				params.add_many(p2)
	
				componentslist.append(comp)
				
				if (comp.name == 'psfwing'):
				
					#psfwing_02pxscale_datatab = readEllipseOutput('PSFtinytim_centered_resc_linscale05px.ell')
					psfwing_02pxscale_datatab = readEllipseOutput('star_02pxscale.ell')
					psfwing_02pxscale_datatab['sma'] = psfwing_02pxscale_datatab['sma'] * Settings.pxlToArcsec
					if equivalentAxisFit:
						psfwing_02pxscale_datatab['sma'] = psfwing_02pxscale_datatab['sma'] * np.sqrt(1 - psfwing_02pxscale_datatab['ellip'])
					#if minorAxisFit:
					#	psfwing_02pxscale_datatab['sma'] = psfwing_02pxscale_datatab['sma'] * (1 - psfwing_02pxscale_datatab['ellip'])
					psfwing_02pxscale_datatab['intens'] = psfwing_02pxscale_datatab['intens'] / Settings.pxlToArcsec**2
					psfwing_02pxscale_datatab['intens'] = psfwing_02pxscale_datatab['intens'] / max(psfwing_02pxscale_datatab['intens'])
					
					#psfwing_logscale_datatab = readEllipseOutput('PSFtinytim_centered_resc_logscale.ell')
					psfwing_logscale_datatab = readEllipseOutput('star_logscale.ell')
					psfwing_logscale_datatab['sma'] = psfwing_logscale_datatab['sma'] * Settings.pxlToArcsec
					if equivalentAxisFit:	
						psfwing_logscale_datatab['sma'] = psfwing_logscale_datatab['sma'] * np.sqrt(1 - psfwing_logscale_datatab['ellip'])
					#if minorAxisFit:
					#	psfwing_logscale_datatab['sma'] = psfwing_logscale_datatab['sma'] * (1 - psfwing_logscale_datatab['ellip'])
					psfwing_logscale_datatab['intens'] = psfwing_logscale_datatab['intens'] / Settings.pxlToArcsec**2
					psfwing_logscale_datatab['intens'] = psfwing_logscale_datatab['intens'] / max(psfwing_logscale_datatab['intens'])
			
	return componentslist, params, psfwing_02pxscale_datatab, psfwing_logscale_datatab

