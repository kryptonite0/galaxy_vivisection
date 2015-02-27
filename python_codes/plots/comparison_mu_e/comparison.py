import sqlite3 as sql3
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.mlab as mlab

from errors.comparison import comparison

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

matplotlib.rcParams.update({'font.size': 22})

#colordict = ['black', 'black', 'blue', 'lightgreen', 'yellow', 'orange']
colordict = ['black', 'blue', 'lightgreen', 'yellow', 'orange']

path_comp_plots = '/Users/gsavorgnan/galaxy_vivisection/results/plots/comparison_mu_e/'
path_paper_images = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/images/'

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

n_grade1 = 27.
n_grade2 = 29.
n_grade3 = 10.
n_tot = 72.
f_grade1 = n_grade1/n_tot
f_grade2 = (n_grade1+n_grade2)/n_tot
f_grade3 = (n_grade1+n_grade2+n_grade3)/n_tot
#print f_grade1,f_grade2,f_grade3

def comparison_mu_e_my():
	
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	query = 'SELECT mu_e_maj_moffat_comb, mu_e_eq_moffat_comb from OneDFitResults;'
	cur.execute(query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	mu_e_maj = data[0].astype(np.float)
	mu_e_eq = data[1].astype(np.float)
	
	fig, ax1 = plt.subplots()
	
	#diff1 = (+mu_e_maj - mu_e_eq)*0.5
	#diff2 = (-mu_e_maj + mu_e_eq)*0.5
	#diff = np.concatenate((diff1, diff2), axis=0)
	diff = mu_e_maj-mu_e_eq
	
	average = np.average(diff)
	sigma = np.std(diff)
	x = np.asarray(np.arange(-3,3,0.01))
	y_gauss = 20*mlab.normpdf( x, average, sigma)
	#ax1.plot(x, y_gauss, 'k--', linewidth=3)
	print 'sigma mu_e = ', sigma
	
	abs_diff = np.sort(abs(diff))
	#print abs_diff
	sigma_grade1 = abs_diff[int(round(len(abs_diff)*f_grade1))-1]
	sigma_grade2 = abs_diff[int(round(len(abs_diff)*f_grade2))-1]
	sigma_grade3 = abs_diff[int(round(len(abs_diff)*f_grade3))-1]
	print 'sigma grade 1,2,3 = ', sigma_grade1,sigma_grade2,sigma_grade3
	
	n, bins, patches = ax1.hist(diff, 15, facecolor='black')
	ax1.set_xlabel(r'$(\mu_{\rm e,maj}-\mu_{\rm e,eq})~\rm [mag]$', labelpad=20)
	ax1.set_ylabel(r'$N$', labelpad=20, rotation=0.1)
	ax1.axis([-1.99,1.99,0.1,22.9])
	plt.subplots_adjust(left=0.15,bottom=0.20)
	
	ax1.errorbar((average, average, average), (17, 19, 21), xerr=([sigma_grade3,sigma_grade2,sigma_grade1],[0,0,0]), color = 'k', linewidth=3, linestyle='None')
	ax1.text(average+0.1, 20.5, r'$\sigma_1 = 0.17~\rm mag$', fontsize=25)
	ax1.text(average+0.1, 18.5, r'$\sigma_2 = 0.56~\rm mag$', fontsize=25)
	ax1.text(average+0.1, 16.5, r'$\sigma_3 = 1.14~\rm mag$', fontsize=25)
	
	#plt.show()
	plt.savefig(path_comp_plots + 'comparison_my_mu_e.eps', format='eps', dpi=1000)	

def comparison_mu_e():
	
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
			measurementslist_mu_e = comparison.getMeasurements_mu_e(cur, gal_id[0])
			#print measurementslist_mu_e
			average_mu_e, stdev_mu_e, num_measurements_mu_e = comparison.computeStatistics(measurementslist_mu_e)
			
			if num_measurements_mu_e>1:
				diff_my = measurementslist_mu_e[0] - average_mu_e
				diff_my_list.append(diff_my)
				#diff_my = measurementslist_mu_e[1] - average_mu_e
				#diff_my_list.append(diff_my)
			
			for i in range(0,len(measurementslist_mu_e)):
				if num_measurements_mu_e>1:
					diff = measurementslist_mu_e[i] - average_mu_e
					diff_list.append(diff)
					ax2.scatter(diff, average_mu_e, c=colordict[i], s=40)
	print 'number of black points = ', len(diff_my_list)					
	diff_array = np.asarray([diff_list])
	diff_array = diff_array[diff_array>-99]
	sigma = np.std(diff_array)
	#print 'sigma mag_sph = ', sigma
	diff_my_array = np.asarray([diff_my_list])
	diff_my_array = diff_my_array[diff_my_array>-99]
	sigma_my = np.std(diff_my_array)
	#print 'my sigma mag_sph = ', sigma_my
	n, bins, patches = ax1.hist(diff_array, 20, facecolor='white')				
	n_my, bins_my, patches_my = ax1.hist(diff_my_array, 10, facecolor='black')				
	ax1.axis([3.99,-3.99,0.1,99])
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
	#print 'sigma grade 1,2,3 = ', sigma_grade1,sigma_grade2,sigma_grade3
	
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

	ax1.errorbar((0, 0, 0), (50, 67.5, 85), xerr=([msigma_grade3,msigma_grade2,msigma_grade1],[psigma_grade3,psigma_grade2,psigma_grade1]), color = 'k', linewidth=3, linestyle='None')
	ax1.text(-msigma_grade1-0.1, 80,   r'$\sigma_1$', fontsize=25)
	ax1.text(-msigma_grade2-0.1, 62.5, r'$\sigma_2$', fontsize=25)
	ax1.text(-msigma_grade3-0.1, 45,   r'$\sigma_3$', fontsize=25)
	
	ax2.axis([3.99,-3.99,20.99,12.01])
	ax2.set_xlabel(r'$\mu_{\rm e,sph} - <\mu_{\rm e,sph}>~{\rm [mag~arcsec^{-2}]}$', labelpad=10)
	ax2.set_ylabel(r'$<\mu_{\rm e,sph}>~{\rm [mag~arcsec^{-2}]}$', labelpad=13)
	
	plt.subplots_adjust(left=0.11,bottom=0.14)
	fig.subplots_adjust(hspace=0)
	#plt.show()
	plt.savefig(path_comp_plots + 'comparison_all_mu_e.eps', format='eps', dpi=1000)
	plt.savefig(path_paper_images + 'comparison_all_mu_e.eps', format='eps', dpi=1000)

	f1 = open(path_comp_plots + 'sigma_mu_e.txt', 'w')
	f1.write('# sign     sigma1    sigma2    sigma3 \n')
	f1.write('+  ' + str(psigma_grade1) + '   ' + str(psigma_grade2) + '   ' + str(psigma_grade3))
	f1.write('\n-  ' + str(msigma_grade1) + '   ' + str(msigma_grade2) + '   ' + str(msigma_grade3))
	f1.close()
	
	
	
def main():
	#comparison_mu_e_my()
	comparison_mu_e()
	
main()		
	
	
