import sqlite3 as sql3
import os
import numpy as np

from conversions import convert 

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
ancillaryfilename = '/Users/gsavorgnan/galaxy_vivisection/data/ancillary/SMBHgalaxysample2013.top'
literaturefilename = '/Users/gsavorgnan/galaxy_vivisection/data/ancillary/literature_decompositions.top'
mosaicsfilename = '/Users/gsavorgnan/galaxy_vivisection/data/galaxies/mosaics.list'

def replaceNull(entry):
	if entry == '""':
		entry = None
	return entry	

def replaceBoolean(entry):
	if entry == '""':
		entry = None
	elif entry == 'True' or entry == 'true':
		entry = 1
	elif entry == 'False' or entry == 'false':			
		entry = 0
	else:
		print 'One value classified as boolean was not of types: true, false or null.'	
	return entry
	
def lookForVelmapAtlas3d(gal_id, fastrotator, slowrotator):
	velmap = None
	status = 0
	if fastrotator or slowrotator:
		filename = '/Users/gsavorgnan/galaxy_vivisection/data/kinematics/' + gal_id + '_velmap_atlas3d.eps'
		if os.path.isfile(filename):
			velmap = open(filename).read()
			status = 1
	return status, velmap		
				
def lookForVelmapSLUGG(gal_id, fastrotator, slowrotator):
	velmap = None
	status = 0
	if fastrotator or slowrotator:
		filename = '/Users/gsavorgnan/galaxy_vivisection/data/kinematics/' + gal_id + '_velmap_slugg.eps'
		if os.path.isfile(filename):
			velmap = open(filename).read()
			status = 1
	return status, velmap		
				
def ancillary():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(ancillaryfilename)

	cur.execute('''CREATE TABLE Ancillary
		(gal_id text, gal_name text, IRAC_image integer, fit1D_done integer, fit2D_done integer, 
		ELLIPTICAL_my integer, edgeon integer, source text, RA real, DEC real, 
		z real, morphtype text, simplemorphtype text, bar integer, distance real, core integer, 
		core_inferred_from_sigma integer, core_my integer, sigma real, mass_BH real, perr_mass_BH real, merr_mass_BH real,
		size_core real, ref_size_core text, KMAG_sph real, KMAG_tot real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_name = line.split()[0]	
			gal_name = replaceNull(gal_name)
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			IRAC_image = line.split()[2]
			IRAC_image = replaceBoolean(IRAC_image)	
			fit1D_done = line.split()[3]
			fit1D_done = replaceBoolean(fit1D_done)	
			fit2D_done = line.split()[4]
			fit2D_done = replaceBoolean(fit2D_done)	
			ELLIPTICAL_my = line.split()[8]
			ELLIPTICAL_my = replaceBoolean(ELLIPTICAL_my)
			core_my = line.split()[9]
			core_my = replaceBoolean(core_my)	 
			edgeon = line.split()[10]
			edgeon = replaceBoolean(edgeon)	 
			source = line.split()[16]
			source = replaceNull(source)
			RA = line.split()[17]
			RA = replaceNull(RA)
			DEC = line.split()[18]
			DEC = replaceNull(DEC)
			z = line.split()[19]
			z = replaceNull(z)
			morphtype = line.split()[20]
			morphtype = replaceNull(morphtype)
			simplemorphtype = line.split()[21]
			simplemorphtype = replaceNull(simplemorphtype)
			bar = line.split()[22]
			bar = replaceBoolean(bar)
			distance = line.split()[23]
			distance = replaceNull(distance)
			core = line.split()[24]
			core = replaceBoolean(core)	
			core_inferred_from_sigma = line.split()[25]
			core_inferred_from_sigma = replaceBoolean(core_inferred_from_sigma)
			size_core = line.split()[26]
			size_core = replaceNull(size_core)
			ref_size_core = line.split()[27]
			ref_size_core = replaceNull(ref_size_core)
			sigma = line.split()[28]
			sigma = replaceNull(sigma)
			mass_BH = line.split()[30]
			mass_BH = replaceNull(mass_BH)
			perr_mass_BH = line.split()[31]
			perr_mass_BH = replaceNull(perr_mass_BH)
			merr_mass_BH = line.split()[32]
			merr_mass_BH = replaceNull(merr_mass_BH)
			KMAG_tot = line.split()[38]
			KMAG_tot = replaceNull(KMAG_tot)
			KMAG_sph = line.split()[47]
			KMAG_sph = replaceNull(KMAG_sph)
			
			collection = [gal_id, gal_name, IRAC_image, fit1D_done, fit2D_done, ELLIPTICAL_my, edgeon, source, RA, DEC, z, morphtype, simplemorphtype, bar, distance, 
				core, core_inferred_from_sigma, core_my, sigma, mass_BH, perr_mass_BH, merr_mass_BH, size_core, ref_size_core, KMAG_sph, KMAG_tot]
			cur.execute('INSERT INTO Ancillary VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', collection)

	connection.commit()
	cur.close()
	connection.close()
	data.close()
	
def kinematics():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(ancillaryfilename)

	cur.execute('''CREATE TABLE Kinematics
		(gal_id text, fastrotator integer, slowrotator integer, velmap_atlas3d integer, 
		velmap_atlas3d_file blob, velmap_slugg integer, velmap_slugg_file blob)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			fastrotator = line.split()[13]
			fastrotator = replaceBoolean(fastrotator)
			slowrotator = line.split()[14]
			slowrotator = replaceBoolean(slowrotator)	
			velmap_atlas3d, velmap_atlas3d_file = lookForVelmapAtlas3d(gal_id, fastrotator, slowrotator)
			velmap_slugg, velmap_slugg_file = lookForVelmapSLUGG(gal_id, fastrotator, slowrotator)
			
			collection = [gal_id, fastrotator, slowrotator, velmap_atlas3d, velmap_atlas3d_file, velmap_slugg, velmap_slugg_file]
			cur.execute('INSERT INTO Kinematics VALUES (?,?,?,?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
	
def literatureGrahamDriver2007():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsGrahamDriver2007
		(gal_id text, n real, r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[6]
			n = replaceNull(n)
			mag_sph = None
			r_e = line.split()[18]
			r_e = replaceNull(r_e)
			
			collection = [gal_id, n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsGrahamDriver2007 VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureGrahamDriver2007PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsGrahamDriver2007PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[6]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = None
			r_e = line.split()[18]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			
			collection = [gal_id, log_n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsGrahamDriver2007PhysicalUnits VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureLaurikainenetal2010():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsLaurikainenetal2010 
		(gal_id text, n real, r_e real, mag_sph real, axis_ratio real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[7]
			n = replaceNull(n)
			mag_sph = line.split()[13]
			mag_sph = replaceNull(mag_sph)
			r_e = line.split()[19]
			r_e = replaceNull(r_e)
			axis_ratio = line.split()[26]
			axis_ratio = replaceNull(axis_ratio)
			
			collection = [gal_id, n, r_e, mag_sph, axis_ratio]		
			cur.execute('INSERT INTO LiteratureDecompositionsLaurikainenetal2010 VALUES (?,?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureLaurikainenetal2010PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsLaurikainenetal2010PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real, axis_ratio real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[7]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = line.split()[13]
			mag_sph = replaceNull(mag_sph)
			if mag_sph is not None:
				mag_sph = convert.apparentToAbsoluteMagnitude(float(mag_sph),dist_Mpc)
			r_e = line.split()[19]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			axis_ratio = line.split()[26]
			axis_ratio = replaceNull(axis_ratio)
			
			collection = [gal_id, log_n, r_e, mag_sph, axis_ratio]		
			cur.execute('INSERT INTO LiteratureDecompositionsLaurikainenetal2010PhysicalUnits VALUES (?,?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureSanietal2011():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsSanietal2011
		(gal_id text, n real, r_e real, mag_sph real, axis_ratio real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[8]
			n = replaceNull(n)
			mag_sph = line.split()[14]
			mag_sph = replaceNull(mag_sph)
			r_e = line.split()[20]
			r_e = replaceNull(r_e)
			axis_ratio = line.split()[25]
			axis_ratio = replaceNull(axis_ratio)
			
			collection = [gal_id, n, r_e, mag_sph, axis_ratio]		
			cur.execute('INSERT INTO LiteratureDecompositionsSanietal2011 VALUES (?,?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureSanietal2011PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsSanietal2011PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real, axis_ratio real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[8]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = line.split()[14]
			mag_sph = replaceNull(mag_sph)
			if mag_sph is not None:
				mag_sph = convert.apparentToAbsoluteMagnitude(float(mag_sph),dist_Mpc)
			r_e = line.split()[20]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			axis_ratio = line.split()[25]
			axis_ratio = replaceNull(axis_ratio)	
			
			collection = [gal_id, log_n, r_e, mag_sph, axis_ratio]		
			cur.execute('INSERT INTO LiteratureDecompositionsSanietal2011PhysicalUnits VALUES (?,?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureVikaetal2012():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsVikaetal2012
		(gal_id text, n real, r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[9]
			n = replaceNull(n)
			mag_sph = line.split()[15]
			mag_sph = replaceNull(mag_sph)
			r_e = line.split()[21]
			r_e = replaceNull(r_e)
			
			collection = [gal_id, n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsVikaetal2012 VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureVikaetal2012PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsVikaetal2012PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[9]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = line.split()[15]
			mag_sph = replaceNull(mag_sph)
			if mag_sph is not None:
				mag_sph = convert.apparentToAbsoluteMagnitude(float(mag_sph),dist_Mpc)
			r_e = line.split()[21]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			
			collection = [gal_id, log_n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsVikaetal2012PhysicalUnits VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureBeifiorietal2012():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsBeifiorietal2012
		(gal_id text, n real, r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[10]
			n = replaceNull(n)
			mag_sph = line.split()[16]
			mag_sph = replaceNull(mag_sph)
			r_e = line.split()[22]
			r_e = replaceNull(r_e)
			
			collection = [gal_id, n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsBeifiorietal2012 VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureBeifiorietal2012PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsBeifiorietal2012PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[10]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = line.split()[16]
			mag_sph = replaceNull(mag_sph)
			if mag_sph is not None:
				mag_sph = convert.apparentToAbsoluteMagnitude(float(mag_sph),dist_Mpc)
			r_e = line.split()[22]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			
			collection = [gal_id, log_n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsBeifiorietal2012PhysicalUnits VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureRuslietal2013():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsRuslietal2013
		(gal_id text, n real, r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[11]
			n = replaceNull(n)
			mag_sph = None
			r_e = line.split()[23]
			r_e = replaceNull(r_e)
			
			collection = [gal_id, n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsRuslietal2013 VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureRuslietal2013PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsRuslietal2013PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[11]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = None
			r_e = line.split()[23]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			
			collection = [gal_id, log_n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsRuslietal2013PhysicalUnits VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
def literatureLaskeretal2014():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsLaskeretal2014
		(gal_id text, n real, r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			n = line.split()[12]
			n = replaceNull(n)
			mag_sph = line.split()[17]
			mag_sph = replaceNull(mag_sph)
			r_e = line.split()[24]
			r_e = replaceNull(r_e)
			
			collection = [gal_id, n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsLaskeretal2014 VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()

def literatureLaskeretal2014PhysicalUnits():
	connection = sql3.connect(dbname)
	cur = connection.cursor()

	data = open(literaturefilename)

	cur.execute('''CREATE TABLE LiteratureDecompositionsLaskeretal2014PhysicalUnits
		(gal_id text, log_n real, log_r_e real, mag_sph real)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[1]	
			gal_id = replaceNull(gal_id)
			#get distance in mpc
			if gal_id is not None:
				cur.execute('SELECT distance FROM Ancillary WHERE gal_id = "' + gal_id + '";')
				dist_Mpc = cur.fetchall()[0][0]

			n = line.split()[12]
			n = replaceNull(n)
			log_n = None
			if n is not None:
				log_n = np.log10(float(n))
			mag_sph = line.split()[17]
			mag_sph = replaceNull(mag_sph)
			if mag_sph is not None:
				mag_sph = convert.apparentToAbsoluteMagnitude(float(mag_sph),dist_Mpc)
			r_e = line.split()[24]
			r_e = replaceNull(r_e)
			if r_e is not None:
				r_e = np.log10(convert.arcsecToKpc(float(r_e),dist_Mpc))
			
			collection = [gal_id, log_n, r_e, mag_sph]		
			cur.execute('INSERT INTO LiteratureDecompositionsLaskeretal2014PhysicalUnits VALUES (?,?,?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
	
def spitzerdata():
	connection = sql3.connect(dbname)
	cur = connection.cursor()
	
	data = open(mosaicsfilename) 

	cur.execute('''CREATE TABLE SpitzerData
		(gal_id text, fits text)''')
	
	for line in data:
		if line.split()[0] != '#':
			gal_id = line.split()[0]	
			gal_id = replaceNull(gal_id)
			fitsname = line.split()[1]
			fitsname = replaceNull(fitsname)
					
			collection = [gal_id, fitsname]		
			cur.execute('INSERT INTO SpitzerData VALUES (?,?)', collection)

	connection.commit()

	cur.close()
	connection.close()
	data.close()
		
		

def main():
	try:
		os.remove(dbname)
	except:
		None
	
	ancillary()
	print 'Table Ancillary has been created.'
	kinematics()
	print 'Table Kinematics has been created.'
	literatureGrahamDriver2007()
	literatureGrahamDriver2007PhysicalUnits()
	print 'Table LiteratureDecompositionsGrahamDriver2007 has been created.'
	literatureLaurikainenetal2010()
	literatureLaurikainenetal2010PhysicalUnits()
	print 'Table LiteratureDecompositionsLaurikainenetal2010 has been created.'
	literatureSanietal2011()
	literatureSanietal2011PhysicalUnits()
	print 'Table LiteratureDecompositionsSanietal2011 has been created.'
	literatureVikaetal2012()
	literatureVikaetal2012PhysicalUnits()
	print 'Table LiteratureDecompositionsVikaetal2012 has been created.'
	literatureBeifiorietal2012()
	literatureBeifiorietal2012PhysicalUnits()
	print 'Table LiteratureDecompositionsBeifiorietal2012 has been created.'
	literatureRuslietal2013()
	literatureRuslietal2013PhysicalUnits()
	print 'Table LiteratureDecompositionsRuslietal2013 has been created.'
	literatureLaskeretal2014()
	literatureLaskeretal2014PhysicalUnits()
	print 'Table LiteratureDecompositionsLaskeretal2014 has been created.'
	#spitzerdata()
	#print 'Table SpitzerData has been created.'
	
	
	
main()	
