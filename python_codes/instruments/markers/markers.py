import numpy as np
import matplotlib
import matplotlib.pyplot as plt

scatter_kwargs = {"zorder":100}

def spiral(ax, color, x0, y0, size):
	t1 = np.arange(1,4,0.01)
	t2 = np.arange(1,4,0.01)
	x1 = size*t1*np.cos(t1)
	y1 = size*t1*np.sin(t1)
	x2 = -size*t2*np.cos(t2)
	y2 = -size*t2*np.sin(t2)
	
	t = np.arange(0,2*np.pi,0.01)
	x = size*1.1*np.cos(t)
	y = size*1*np.sin(t)
	
	ax.plot(10**(x0+x1),10**(y0+y1), color=color, lw=2, **scatter_kwargs)
	ax.plot(10**(x0+x2),10**(y0+y2), color=color, lw=2, **scatter_kwargs)
	ax.plot(10**(x0+x),10**(y0+y), color=color, lw=2, **scatter_kwargs)
	ax.fill_between(10**(x0+x), 10**(y0+0), 10**(y0+y), edgecolor='', facecolor=color, **scatter_kwargs)
	
	
def lenticular(ax, color, x0, y0, size):

	t = np.arange(0,2*np.pi,0.01)
	x = size*1*np.cos(t)
	y = size*1.5*np.sin(t)	
	
	ax.fill_between(10**(x0+x), 10**(y0+0), 10**(y0+y), edgecolor='', facecolor=color, **scatter_kwargs)
	ax.plot([10**(x0-2*size),10**(x0+2*size)],[10**(y0),10**(y0)], color=color, lw=3, **scatter_kwargs)


def elliptical(ax, color, x0, y0, size):

	t = np.arange(0,2*np.pi,0.01)
	x = size*1*np.cos(t)
	y = size*1.5*np.sin(t)	
	
	ax.fill_between(10**(x0+x), 10**(y0+0), 10**(y0+y), edgecolor='', facecolor=color, **scatter_kwargs)


