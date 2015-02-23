import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pyfits

from scipy.misc import imread
import matplotlib.cbook as cbook

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 15})

pxlToArcsec = 1.22
fwhm_gauss = 2.03/pxlToArcsec
alpha_moff = 2.38/pxlToArcsec
beta_moff = 4.39

pathinput = '/Users/gsavorgnan/galaxy_vivisection/python_codes/paper_figures/datapaper/psf/'
pathoutput = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/images/'

def getRadialProfile(fitsFileName, xc, yc):
	
	f = pyfits.open(pathinput + fitsFileName)  # open a FITS file
	scidata = f[0].data  # assume the first extension is an image
	radius = [0.0]*scidata.size
	intens = [0.0]*scidata.size
	
	for i in range(0,scidata.shape[0]):
		for j in range(0,scidata.shape[1]):
			R = ((i-xc)**2+(j-yc)**2)**0.5
			I = scidata[i,j]
			radius[i+j*scidata.shape[0]] = R
			intens[i+j*scidata.shape[0]] = I
	
	return radius, intens
	
def plotRadialProfile(fig, ax, fitsFileName, imageFileName, title, ypos, xc, yc, rebin):
	
	radius, intens = getRadialProfile(fitsFileName, xc, yc)
	
	radius = np.asarray(radius)
	intens = np.asarray(intens)
	
	radius_red = radius[radius<6*rebin]
	intens_red = intens[radius<6*rebin]
	
	radius_pixel = radius_red / rebin
	intens_norm = intens_red/max(intens_red)
			
	#fig, ax1 = plt.subplots()
	ax.scatter(radius_pixel, intens_norm, c='black', s=40, zorder=1)
        ax.axis([-0.9,(6.9),-0.1,1.1])
        #ax.set_xlabel(r'$R \rm~[pixel]$', labelpad=10)
        ax.set_ylabel(r'norm. flux', labelpad=10)
	#plt.title(title)
	ax.text(3.75, 0.25, title, fontsize=20)
	#plt.subplots_adjust(left=0.20,bottom=0.20)
	
	xxx = np.arange(0,10,0.05)
	A = 1/np.exp(-(min(radius_pixel))**2/(2*(fwhm_gauss/2.355)**2))
	yyy_gauss = A*np.exp(-xxx**2/(2*(fwhm_gauss/2.355)**2))
	B = 1/(1+(min(radius_pixel)/alpha_moff)**2)**(-beta_moff)
	yyy_moff = B*(1+(xxx/alpha_moff)**2)**(-beta_moff)
	ax.plot(xxx, yyy_moff, c='blue', linewidth=3)
	#plt.plot(xxx, yyy_gauss, c='red', linewidth=3, ls='--')
	
	axplot = fig.add_axes([0.72,ypos,0.17,0.17])
	axplot.set_xticks([])
	axplot.set_yticks([])
	img = imread(pathinput + imageFileName)
	axplot.imshow(img)
	
	#plt.show()
	#plt.savefig(pathoutput+fitsFileName.replace('fits','eps'), format='eps', dpi=1000)
	return fig, ax
	
def main():
	
	fig, (ax1, ax2) = plt.subplots(2, sharex=True)
	fig, ax1 = plotRadialProfile(fig, ax1, 'PRF_IRACch1_129_129_cut.fits', 'IRACPRF.jpeg', r'IRAC PRF', 0.7, 49.67, 49.54, 5)
	fig, ax2 = plotRadialProfile(fig, ax2, 'star1_cut.fits', 'star1.jpeg', r'Real star', 0.35, 9.86, 9.69, 1)
	ax2.set_xlabel(r'$R \rm~[pixel]$', labelpad=10)
	fig.subplots_adjust(left=0.50,bottom=0.20)
	fig.subplots_adjust(hspace=0)
	
	#plt.show()
	plt.savefig(pathoutput+'psf.eps', format='eps', dpi=1000)
	
main()
