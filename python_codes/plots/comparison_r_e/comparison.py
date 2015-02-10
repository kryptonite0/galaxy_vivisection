import sqlite3 as sql3
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.mlab as mlab

from errors.comparison import comparison

# exclude R_eq of rusli, the others are R_maj

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 22})

colordict = ['black', 'red', 'blue', 'lightgreen', 'yellow', 'gray', 'orange']

path_comp_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_r_e/'
path_paper_images = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/images/'

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

#workslist = ['GrahamDriver2007', 'Laurikainenetal2010', 'Sanietal2011', 'Vikaetal2012',
#	'Beifiorietal2012', 'Laskeretal2014']

n_grade1 = 27.
n_grade2 = 29.
n_grade3 = 10.
n_tot = 72.
f_grade1 = n_grade1/n_tot
f_grade2 = (n_grade1+n_grade2)/n_tot
f_grade3 = (n_grade1+n_grade2+n_grade3)/n_tot
#print f_grade1,f_grade2,f_grade3

def comparison_log_r_e():
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	query = 'select gal_id from OneDFitResults;'
	cur.execute(query)
	galaxieslist = cur.fetchall()
	
	fig, (ax1, ax2) = plt.subplots(2, sharex=True)

	diff_list = []
	diff_my_list = []
	for gal_id in galaxieslist:
		if gal_id[0] is not None:
			measurementslist_log_r_e = comparison.getMeasurements_log_r_e(cur, gal_id[0])
			average_log_r_e, stdev_log_r_e, num_measurements_log_r_e = comparison.computeStatistics(measurementslist_log_r_e)
			
			if num_measurements_log_r_e>1:
				diff_my = measurementslist_log_r_e[0] - average_log_r_e
				diff_my_list.append(diff_my)
			
			for i in range(0,len(measurementslist_log_r_e)):
				if num_measurements_log_r_e>1:
					diff = measurementslist_log_r_e[i] - average_log_r_e
					diff_list.append(diff)
					ax2.scatter(diff, average_log_r_e, c=colordict[i], s=40)
	print 'number of black points = ', len(diff_my_list)					
	diff_array = np.asarray([diff_list])
	diff_array = diff_array[diff_array>-99]
	sigma = np.std(diff_array)
	#print 'sigma log_r_e = ', sigma
	diff_my_array = np.asarray([diff_my_list])
	diff_my_array = diff_my_array[diff_my_array>-99]
	sigma_my = np.std(diff_my_array)
	#print 'my sigma log_r_e = ', sigma_my
	n, bins, patches = ax1.hist(diff_array, 15, facecolor='white')				
	n_my, bins_my, patches_my = ax1.hist(diff_my_array, 10, facecolor='black')				
	ax1.axis([-1.29,1.29,0.1,99])
	x = np.asarray(np.arange(-3,3,0.01))
	y_gauss = 27*mlab.normpdf( x, 0, sigma)
	y_gauss_my = 8*mlab.normpdf( x, 0, sigma_my)
	#ax1.plot(x, y_gauss, color='black', linestyle='-', linewidth=3)
	#ax1.plot(x, y_gauss_my, 'k--', linewidth=3)
	ax1.set_ylabel(r'$N$', labelpad=20, rotation=0.1)
		
	#abs_diff_my_array = np.sort(abs(diff_my_array))
	#sigma_grade1 = abs_diff_my_array[int(len(abs_diff_my_array)*f_grade1)]
	#sigma_grade2 = abs_diff_my_array[int(len(abs_diff_my_array)*f_grade2)]
	#sigma_grade3 = abs_diff_my_array[int(len(abs_diff_my_array)*f_grade3)]
	
	pdiff_my_array = diff_my_array[diff_my_array>=0]
	mdiff_my_array = abs(diff_my_array[diff_my_array<0])
	pdiff_my_array_sorted = np.sort(pdiff_my_array)
	mdiff_my_array_sorted = np.sort(mdiff_my_array)
	#print pdiff_my_array_sorted,mdiff_my_array_sorted
	
	psigma_grade1 = pdiff_my_array_sorted[int(round(len(pdiff_my_array_sorted)*f_grade1))-1]
	psigma_grade2 = pdiff_my_array_sorted[int(round(len(pdiff_my_array_sorted)*f_grade2))-1]
	psigma_grade3 = pdiff_my_array_sorted[int(round(len(pdiff_my_array_sorted)*f_grade3))-1]	
	print '+sigma grade 1,2,3 = ', psigma_grade1,psigma_grade2,psigma_grade3
	
	msigma_grade1 = mdiff_my_array_sorted[int(round(len(mdiff_my_array_sorted)*f_grade1))-1]
	msigma_grade2 = mdiff_my_array_sorted[int(round(len(mdiff_my_array_sorted)*f_grade2))-1]
	msigma_grade3 = mdiff_my_array_sorted[int(round(len(mdiff_my_array_sorted)*f_grade3))-1]
	print '-sigma grade 1,2,3 = ', msigma_grade1,msigma_grade2,msigma_grade3

	ax1.errorbar((0, 0, 0), (60, 72.5, 85), xerr=([msigma_grade3,msigma_grade2,msigma_grade1],[psigma_grade3,psigma_grade2,psigma_grade1]), color = 'k', linewidth=3, linestyle='None')
	ax1.text(psigma_grade1+0.06, 80,   r'$\sigma_1$', fontsize=25)
	ax1.text(psigma_grade2+0.06, 67.5, r'$\sigma_2$', fontsize=25)
	ax1.text(psigma_grade3+0.06, 55,   r'$\sigma_3$', fontsize=25)
	
	ax2.axis([-1.29,1.29,np.log10(1.01),np.log10(500)])
	ax2.set_xlabel(r'$\log(R_{\rm e,sph}^{\rm maj} ~{\rm [arcsec]}) - <\log(R_{\rm e,sph}^{\rm maj}~{\rm [arcsec]})>$', labelpad=20)
	ax2.set_ylabel(r'$<\log(R_{\rm e,sph}^{\rm maj} ~{\rm [arcsec]})>$', labelpad=20)
	
	plt.subplots_adjust(left=0.15,bottom=0.17)
	fig.subplots_adjust(hspace=0)
	#plt.show()
	plt.savefig(path_comp_plots + 'comparison_all_r_e.eps', format='eps', dpi=1000)
	plt.savefig(path_paper_images + 'comparison_all_r_e.eps', format='eps', dpi=1000)
	
	f1 = open(path_comp_plots + 'sigma_log_r_e.txt', 'w')
	f1.write('# sign     sigma1    sigma2    sigma3 \n')
	f1.write('+  ' + str(psigma_grade1) + '   ' + str(psigma_grade2) + '   ' + str(psigma_grade3))
	f1.write('\n-  ' + str(msigma_grade1) + '   ' + str(msigma_grade2) + '   ' + str(msigma_grade3))
	f1.close()
	
	
	
def main():
	comparison_log_r_e()

main()		
	
	
