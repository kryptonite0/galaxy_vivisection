import numpy as np

def arcsecToKpc(arcsec, dist_Mpc):
	dist_Kpc = dist_Mpc * 1000.
	kpc = (arcsec / 3600.) * (np.pi/180.) * dist_Kpc
	return kpc

def apparentToAbsoluteMagnitude(apparent, dist_Mpc):
	absolute = apparent - 5. * np.log10(dist_Mpc) - 25.
	return absolute
