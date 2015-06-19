import numpy as np
import matplotlib
import matplotlib.pyplot as plt

scatter_kwargs = {"zorder":100}
scatter_kwargs2 = {"zorder":200}


def spiral(ax, color, x0, y0, size1, size2):
	t1 = np.arange(1,4,0.01)
	t2 = np.arange(1,4,0.01)
	x1 = size1*t1*np.cos(t1)
	y1 = size2*t1*np.sin(t1)
	x2 = -size1*t2*np.cos(t2)
	y2 = -size2*t2*np.sin(t2)
	
	t = np.arange(0,2*np.pi,0.01)
	x = size1*1.1*np.cos(t)
	y = size2*1*np.sin(t)
	
	ax.plot((x0+x1),(y0+y1), color=color, lw=2, **scatter_kwargs2)
	ax.plot((x0+x2),(y0+y2), color=color, lw=2, **scatter_kwargs2)
	ax.plot((x0+x), (y0+y), color=color, lw=2, **scatter_kwargs2)
	ax.fill_between((x0+x), (y0+0), (y0+y), edgecolor='', facecolor=color, **scatter_kwargs2)
	
	
def lenticular(ax, color, x0, y0, size1, size2):

	t = np.arange(0,2*np.pi,0.01)
	x = size1*1*np.cos(t)
	y = size2*1.5*np.sin(t)	
	
	ax.fill_between((x0+x), (y0+0), (y0+y), edgecolor='', facecolor=color, **scatter_kwargs)
	ax.plot([(x0-2*size1),(x0+2*size1)],[(y0),(y0)], color=color, lw=3, **scatter_kwargs)


def elliptical(ax, color, x0, y0, size1, size2):

	t = np.arange(0,2*np.pi,0.01)
	x = size1*2*np.cos(t)
	y = size2*1.5*np.sin(t)	
	
	ax.fill_between((x0+x), (y0+0), (y0+y), edgecolor='', facecolor=color, **scatter_kwargs)


