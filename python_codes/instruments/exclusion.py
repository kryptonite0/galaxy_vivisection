import numpy as np

def clip(a, threshold):
	if len(a)==1:
		return a
	if len(a)>1:
		a.sort()
		mask = np.ones(a.size, bool)
		for i in range(len(a)-1):
			diff = (a[i+1]-a[i])
			if diff>threshold:
				mask[i] = False
		last_diff = a[-1] - a[-2]
		if last_diff >threshold:
			mask[-1] = False		
		return a[mask]		
