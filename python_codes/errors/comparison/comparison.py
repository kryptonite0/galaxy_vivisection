import sqlite3 as sql3
import os
import numpy as np

from scipy.special import gamma
import b_n

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

workslist_r_e = ['GrahamDriver2007', 'Laurikainenetal2010', 'Sanietal2011', 'Vikaetal2012',
	'Beifiorietal2012', 'Laskeretal2014']
workslist_n = ['GrahamDriver2007', 'Laurikainenetal2010', 'Sanietal2011', 'Vikaetal2012',
	'Beifiorietal2012', 'Ruslietal2013', 'Laskeretal2014']
	
# conversion from 3.6um magnitudes to K-band magnitudes, assumin Salpeter IMF, single stellar population, 13Gyr, solar metallicity
spitzer36ToKbandMag = +0.27	

parameterslist = ['n', 'r_e']

def main():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	try:
		cur.execute('DROP TABLE ErrorsComparison')
		print 'Table ErrorsComparison has been erased.'
	except:
		print 'Table ErrorsComparison does not exist.'
	
	
	create_table_string = 'CREATE TABLE ErrorsComparison (gal_id text, '
	for parameter in parameterslist:
		 create_table_string += ('log_' + parameter + '_average real, ')
		 create_table_string += ('log_' + parameter + '_std real, ') 
		 create_table_string += ('log_' + parameter + '_num_measurements integer, ')
	#create_table_string = create_table_string[:-2]
	create_table_string += 'mag_sph_average real, mag_sph_std real, mag_sph_num_measurements integer, '
	create_table_string += 'mu_e_average real, mu_e_std real, mu_e_num_measurements integer, '
	create_table_string += 'mu_0_average real, mu_0_std real, mu_0_num_measurements integer)'
	#create_table_string += ')'
	cur.execute(create_table_string)
	
	query = 'select gal_id from OneDFitResults;'
	
	cur.execute(query)
	galaxieslist = cur.fetchall()
		
	for gal_id in galaxieslist:
		collection = []
		collection.append(gal_id[0])
		
		measurementslist_log_n = getMeasurements_log_n(cur, gal_id[0])
		average_log_n, stdev_log_n, num_measurements_log_n = computeStatistics(measurementslist_log_n)
		collection.append(average_log_n)
		collection.append(stdev_log_n)
		collection.append(num_measurements_log_n)
		#print gal_id[0], measurementslist_n, average_n, stdev_n, num_measurements_n, ' --- ',
		
		measurementslist_log_r_e = getMeasurements_log_r_e(cur, gal_id[0])
		average_log_r_e, stdev_log_r_e, num_measurements_log_r_e = computeStatistics(measurementslist_log_r_e)
		collection.append(average_log_r_e)
		collection.append(stdev_log_r_e)
		collection.append(num_measurements_log_r_e)
		#print measurementslist_r_e, average_r_e, stdev_r_e, num_measurements_r_e
					
		measurementslist_mag_sph = getMeasurements_mag_sph(cur, gal_id[0])
		average_mag_sph, stdev_mag_sph, num_measurements_mag_sph = computeStatistics(measurementslist_mag_sph)
		collection.append(average_mag_sph)
		collection.append(stdev_mag_sph)
		collection.append(num_measurements_mag_sph)
		#print measurementslist_r_e, average_r_e, stdev_r_e, num_measurements_r_e
		
		measurementslist_mu_e = getMeasurements_mu_e(cur, gal_id[0])
		average_mu_e, stdev_mu_e, num_measurements_mu_e = computeStatistics(measurementslist_mu_e)
		collection.append(average_mu_e)
		collection.append(stdev_mu_e)
		collection.append(num_measurements_mu_e)
					
		measurementslist_mu_0 = getMeasurements_mu_0(cur, gal_id[0])
		average_mu_0, stdev_mu_0, num_measurements_mu_0 = computeStatistics(measurementslist_mu_0)
		collection.append(average_mu_0)
		collection.append(stdev_mu_0)
		collection.append(num_measurements_mu_0)
					
		question_marks = ''
		for i in range(len(collection)-1):
			question_marks += '?,'
		question_marks += '?'	
		
		#print 'INSERT INTO ErrorsComparison VALUES (' + question_marks + ')', collection
		cur.execute('INSERT INTO ErrorsComparison VALUES (' + question_marks + ')', collection)

	connection.commit()
	cur.close()
	connection.close()
	
	print 'Table ErrorsComparison has been created.'

def computeStatistics(lista):
	a = np.asarray(lista, dtype=np.float)
	a = a[np.isfinite(a)]
	num_measurements = len(a)
	average = np.average(a)
	#median = np.median(a)
	stdev = np.std(a)
	#print a, median, stdev, num_measurements
	
	#return median, stdev, num_measurements
	return average, stdev, num_measurements
	
def getMeasurements_log_n(cur, gal_id):	
	
	measurementslist_log_n = []
	
	query_my = 'SELECT n_maj_moffat_comb FROM OneDFitResults WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		log_measurement = np.log10(measurement)
	else:
		log_measurement = np.nan	
	measurementslist_log_n.append(log_measurement)
	
	for work in workslist_n:
		tablename = 'LiteratureDecompositions' + work 
		query_other = 'SELECT n FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
		#print query_other
		cur.execute(query_other)
		data = cur.fetchall()[0]
		measurement = data[0]
		
		if measurement is not None:
			log_measurement = np.log10(measurement)
		else:
			log_measurement = np.nan
		measurementslist_log_n.append(log_measurement)
		
	#print gal_id, parameter, measurementslist
	return measurementslist_log_n
	
def getMeasurements_log_r_e(cur, gal_id):	
	
	measurementslist_log_r_e = []
	
	query_my = 'SELECT r_e_maj_moffat_comb FROM OneDFitResults WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		log_measurement = np.log10(measurement)
	else:
		log_measurement = np.nan	
	measurementslist_log_r_e.append(log_measurement)
	
	for work in workslist_r_e:
		tablename = 'LiteratureDecompositions' + work 
		query_other = 'SELECT r_e FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
		#print query_other
		cur.execute(query_other)
		#print cur.fetchall(), gal_id, parameter
		measurement = cur.fetchall()[0][0]
		if measurement is not None:
			log_measurement = np.log10(measurement)
		else:
			log_measurement = np.nan
		measurementslist_log_r_e.append(log_measurement)
		
	#print gal_id, parameter, measurementslist
	return measurementslist_log_r_e
	
def getMeasurements_mag_sph(cur, gal_id):	
	
	measurementslist_mag_sph = []
	
	query_my = 'SELECT mag_sph_eq_moffat_comb FROM OneDFitResults WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		measurement = measurement + spitzer36ToKbandMag
	else:
		measurement = np.nan	
	measurementslist_mag_sph.append(measurement)
	
	# get mags from laurikainen (K-band)
	tablename = 'LiteratureDecompositionsLaurikainenetal2010'
	query_other = 'SELECT mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		measurement = measurement
	else:
		measurement = np.nan	
	measurementslist_mag_sph.append(measurement)
	
	# get mags from sani (3.6-band)
	tablename = 'LiteratureDecompositionsSanietal2011'
	query_other = 'SELECT mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		measurement = measurement + spitzer36ToKbandMag
	else:
		measurement = np.nan	
	measurementslist_mag_sph.append(measurement)
	
	# get mags from vika (K-band)
	tablename = 'LiteratureDecompositionsVikaetal2012'
	query_other = 'SELECT mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		measurement = measurement
	else:
		measurement = np.nan	
	measurementslist_mag_sph.append(measurement)
	
	# get mags from lasker (K-band)
	tablename = 'LiteratureDecompositionsLaskeretal2014'
	query_other = 'SELECT mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	measurement = cur.fetchall()[0][0]
	if measurement is not None:
		measurement = measurement
	else:
		measurement = np.nan	
	measurementslist_mag_sph.append(measurement)
		
	#print gal_id, parameter, measurementslist
	return measurementslist_mag_sph
	
def getMeasurements_mu_e(cur, gal_id):	
	
	measurementslist_mu_e = []
	
	# get mu_e from my sample
	query_my = 'SELECT mu_e_maj_moffat_comb, mu_e_eq_moffat_comb FROM OneDFitResults WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	mu_e_maj = data[0].astype(np.float)
	if mu_e_maj is not None:
		mu_e_maj = mu_e_maj 
	else:
		mu_e_maj = np.nan	
	measurementslist_mu_e.append(mu_e_maj)
	
	mu_e_eq  = data[1].astype(np.float)
	if mu_e_eq is not None:
		mu_e_eq = mu_e_eq 
	else:
		mu_e_eq = np.nan	
	#measurementslist_mu_e.append(mu_e_eq)
	
	# get mu_e from laurikainen (K-band)
	tablename = 'LiteratureDecompositionsLaurikainenetal2010'
	query_other = 'SELECT n, r_e, mag_sph, axis_ratio FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	n = data[0].astype(np.float)
	r_e = data[1].astype(np.float)
	mag_sph = data[2].astype(np.float)
	axis_ratio = data[3].astype(np.float)
	
	if n is not None and r_e is not None and mag_sph is not None and axis_ratio is not None:
		b = b_n.computeb_n(n) # 1.9992*n-0.3271
		mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
	else:
		mu_e = np.nan
		
	measurementslist_mu_e.append(mu_e)
	
	# get mu_e from sani (3.6-band)
	tablename = 'LiteratureDecompositionsSanietal2011'
	query_other = 'SELECT n, r_e, mag_sph, axis_ratio FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	n = data[0].astype(np.float)
	r_e = data[1].astype(np.float)
	mag_sph = data[2].astype(np.float)
	axis_ratio_sani = data[3].astype(np.float)

	if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
		b = b_n.computeb_n(n) # 1.9992*n-0.3271
		#print b,n
		mu_e = (mag_sph) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
	else:
		mu_e = np.nan
		
	measurementslist_mu_e.append(mu_e)
		
        # get mu_e from vika (K-band)
        tablename = 'LiteratureDecompositionsVikaetal2012'
        query_other = 'SELECT n, r_e, mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
        cur.execute(query_other)
        datalist = cur.fetchall()
        data = np.asarray(datalist).transpose()
       
        n = data[0].astype(np.float)
        r_e = data[1].astype(np.float)
        mag_sph = data[2].astype(np.float)
       
        if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
        	b = b_n.computeb_n(n) # 1.9992*n-0.3271
        	#print b,n
        	mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
        else:
        	mu_e = np.nan
        	
        measurementslist_mu_e.append(mu_e)
       
        # get mu_e from lasker (K-band)
        tablename = 'LiteratureDecompositionsLaskeretal2014'
        query_other = 'SELECT n, r_e, mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
        cur.execute(query_other)
        datalist = cur.fetchall()
        data = np.asarray(datalist).transpose()
       
        n = data[0].astype(np.float)
        r_e = data[1].astype(np.float)
        mag_sph = data[2].astype(np.float)
       
        if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
        	b = b_n.computeb_n(n) # 1.9992*n-0.3271
        	#print b,n
        	mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
        else:
        	mu_e = np.nan
        	
        measurementslist_mu_e.append(mu_e)
       
        return measurementslist_mu_e
	
def getMeasurements_mu_0(cur, gal_id):	
	
	measurementslist_mu_0 = []
	
	# get mu_e from my sample
	query_my = 'SELECT n_maj_moffat_comb, n_eq_moffat_comb, mu_e_maj_moffat_comb, mu_e_eq_moffat_comb FROM OneDFitResults WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	n_maj = data[0].astype(np.float)
	mu_e_maj = data[2].astype(np.float)
	if n_maj is not None and mu_e_maj is not None:
		b_maj = b_n.computeb_n(n_maj)
		mu_0_maj = mu_e_maj - 2.5*np.log10(2*np.pi*n_maj*(np.exp(b_maj)/b_maj**(2*n_maj)) * gamma(2*n_maj))
	else:
		mu_0_maj = np.nan	
	measurementslist_mu_0.append(mu_0_maj)
	
	n_eq = data[1].astype(np.float)
	mu_e_eq = data[3].astype(np.float)
	if n_eq is not None and mu_e_eq is not None:
		b_eq = b_n.computeb_n(n_eq)
		mu_0_eq = mu_e_eq - 2.5*np.log10(2*np.pi*n_eq*(np.exp(b_eq)/b_eq**(2*n_eq)) * gamma(2*n_eq))
	else:
		mu_0_eq = np.nan	
	#measurementslist_mu_0.append(mu_0_eq)
	
	# get mu_e from laurikainen (K-band)
	tablename = 'LiteratureDecompositionsLaurikainenetal2010'
	query_other = 'SELECT n, r_e, mag_sph, axis_ratio FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	n = data[0].astype(np.float)
	r_e = data[1].astype(np.float)
	mag_sph = data[2].astype(np.float)
	axis_ratio = data[3].astype(np.float)
	
	if n is not None and r_e is not None and mag_sph is not None and axis_ratio is not None:
		b = b_n.computeb_n(n) # 1.9992*n-0.3271
		mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))
		mu_0 = mu_e - 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
	else:
		mu_0 = np.nan
		
	measurementslist_mu_0.append(mu_0)
	
	# get mu_e from sani (3.6-band)
	tablename = 'LiteratureDecompositionsSanietal2011'
	query_other = 'SELECT n, r_e, mag_sph, axis_ratio FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_other)
	datalist = cur.fetchall()
	data = np.asarray(datalist).transpose()
	
	n = data[0].astype(np.float)
	r_e = data[1].astype(np.float)
	mag_sph = data[2].astype(np.float)
	axis_ratio_sani = data[3].astype(np.float)

	if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
		b = b_n.computeb_n(n) # 1.9992*n-0.3271
		mu_e = (mag_sph) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
		mu_0 = mu_e - 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
	else:
		mu_0 = np.nan
		
	measurementslist_mu_0.append(mu_0)
		
        # get mu_e from vika (K-band)
        tablename = 'LiteratureDecompositionsVikaetal2012'
        query_other = 'SELECT n, r_e, mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
        cur.execute(query_other)
        datalist = cur.fetchall()
        data = np.asarray(datalist).transpose()
       
        n = data[0].astype(np.float)
        r_e = data[1].astype(np.float)
        mag_sph = data[2].astype(np.float)
       
        if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
        	b = b_n.computeb_n(n) # 1.9992*n-0.3271
        	mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
		mu_0 = mu_e - 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
        else:
        	mu_0 = np.nan
        	
        measurementslist_mu_0.append(mu_0)
       
        # get mu_e from lasker (K-band)
        tablename = 'LiteratureDecompositionsLaskeretal2014'
        query_other = 'SELECT n, r_e, mag_sph FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
        cur.execute(query_other)
        datalist = cur.fetchall()
        data = np.asarray(datalist).transpose()
       
        n = data[0].astype(np.float)
        r_e = data[1].astype(np.float)
        mag_sph = data[2].astype(np.float)
       
        if n is not None and r_e is not None and mag_sph is not None and axis_ratio_sani is not None:
        	b = b_n.computeb_n(n) # 1.9992*n-0.3271
        	mu_e = (mag_sph - spitzer36ToKbandMag) + 5*np.log10(r_e*axis_ratio_sani**0.5) + 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
		mu_0 = mu_e - 2.5*np.log10(2*np.pi*n*(np.exp(b)/b**(2*n)) * gamma(2*n))		
        else:
        	mu_0 = np.nan
        	
        measurementslist_mu_0.append(mu_0)
       
	return measurementslist_mu_0
		

main()
	
##comparison_r_e()	
		
	
		
