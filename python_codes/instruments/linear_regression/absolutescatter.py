import numpy as np

def get_absscatter(x, y, afit, bfit):
	N = len(x)
	Ndof = N-2
	absscat = 0.0
	for i in range(0,N):
		absscat = absscat + (y[i] - afit - bfit*x[i])**2
	absscat = (absscat/Ndof)**0.5 
	
	return absscat
