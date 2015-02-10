import sqlite3 as sql3
import os
import numpy as np

from instruments import exclusion

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

workslist = ['GrahamDriver2007', 'Laurikainenetal2010', 'Sanietal2011', 'Vikaetal2012',
	'Beifiorietal2012', 'Ruslietal2013', 'Laskeretal2014']
	
parameterslist = ['log_n', 'log_r_e', 'mag_sph']
threshold_dict = {'log_n' : np.log10(1.5), 'log_r_e' : np.log10(1.66), 'mag_sph' : np.log10(2.)}	

def main():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	try:
		cur.execute('DROP TABLE ErrorsComparisonSigmaClip')
		print 'Table ErrorsComparisonSigmaClip has been erased.'
	except:
		print 'Table ErrorsComparisonSigmaClip does not exist.'
	
	
	create_table_string = 'CREATE TABLE ErrorsComparisonSigmaClip (gal_id text, '
	for parameter in parameterslist:
		 create_table_string += (parameter + '_median real, ')
		 create_table_string += (parameter + '_std real, ') 
		 create_table_string += (parameter + '_num_measurements integer, ')
	create_table_string = create_table_string[:-2]
	create_table_string += ')'
	cur.execute(create_table_string)
	
	query = 'select gal_id from OneDFitResultsPhysicalUnits;'
	cur.execute(query)
	galaxieslist = cur.fetchall()
		
	for gal_id in galaxieslist:
		collection = []
		collection.append(gal_id[0])
		for parameter in parameterslist:
			measurementslist = getMeasurements(cur, parameter, gal_id[0])
			#print measurementslist
			threshold = threshold_dict[parameter]
			median, stdev, num_measurements = computeStatistics(measurementslist, threshold)
			
			collection.append(median)
			collection.append(stdev)
			collection.append(num_measurements)
			
		question_marks = ''
		for i in range(len(collection)-1):
			question_marks += '?,'#
		question_marks += '?'	
		
		#print 'INSERT INTO ErrorsComparisonSigmaClip VALUES (' + question_marks + ')', collection
		cur.execute('INSERT INTO ErrorsComparisonSigmaClip VALUES (' + question_marks + ')', collection)

	connection.commit()
	cur.close()
	connection.close()
	
	print 'Table ErrorsComparisonSigmaClip has been created.'

def computeStatistics(lista, threshold):
	a = np.asarray(lista, dtype=np.float)
	a = a[np.isfinite(a)]
	a_clipped = exclusion.clip(a, threshold)
	#print a, a_clipped, threshold
	num_measurements = len(a_clipped)
	#average = np.average(a)
	if len(a)>0:
		median = np.median(a_clipped)
		stdev = np.std(a_clipped)
	else:
		median = None
		stdev = None	
	#if len(a) > len(a_clipped):
	#	print len(a), num_measurements
	#print a, a_clipped, threshold, median, stdev, num_measurements
	return median, stdev, num_measurements
	
def getMeasurements(cur, parameter, gal_id):	
	
	measurementslist = []
	
	query_my = 'SELECT ' + parameter + '_eq_moffat_comb FROM OneDFitResultsPhysicalUnits WHERE gal_id = "' + gal_id + '";' 
	cur.execute(query_my)
	measurement = cur.fetchall()[0][0]
	measurementslist.append(measurement)
	
	for work in workslist:
		tablename = 'LiteratureDecompositions' + work + 'PhysicalUnits'
		query_other = 'SELECT ' + parameter + ' FROM ' + tablename + ' WHERE gal_id = "' + gal_id + '";' 
		#print query_other
		cur.execute(query_other)
		#print cur.fetchall(), gal_id, parameter
		measurement = cur.fetchall()[0][0]
		measurementslist.append(measurement)
		
	#print gal_id, parameter, measurementslist
	return measurementslist
	
	
main()

