import sqlite3 as sql3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from scipy import stats
from instruments.linear_regression import bces
from instruments.linear_regression import predband
from instruments.linear_regression import fitexy
from instruments import b_n
from instruments.linear_regression import colorline
from instruments.linear_regression import akaike
import sys

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CenturyGothic']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


matplotlib.rcParams.update({'font.size': 22})

dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'

path_tables = '/Users/gsavorgnan/galaxy_vivisection/results/tables/'

def make_richard_table():
	connection = sql3.connect(dbname)
        cur = connection.cursor()

        cur.execute('''SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.core_inferred_from_sigma, anc.distance,  
		anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, 
		physres.mag_sph_eq_moffat_comb, 
		errV.perr_mag_sph, errV.merr_mag_sph, 
		physres.mag_tot_eq_moffat_comb, 
		col.color, anc.bar,
		res.delta_eq_moffat_comb
		FROM Ancillary AS anc 
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id 
		JOIN OneDFitResults AS res ON anc.gal_id = res.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id 
		JOIN Colors as col ON anc.gal_id = col.gal_id 
		WHERE anc.fit1D_done = 1
                ORDER BY anc.gal_id;''')
                
        datat = cur.fetchall()
        data= np.asarray(datat).transpose()
	
        shape = np.shape(data)
        columns = shape[0]
        rows = shape[1]
        #for i in range(columns):
        #        for j in range(rows):
        #                data[i,j] = putBlankInPlaceOfNone(data[i,j])
        
        gal_name = data[0]
	gal_id = data[0]
	morphtype = data[1]
        core = data[2]
        core[core=='1'] = 'yes'
        core[core=='0'] = 'no'
        core_inferred_from_sigma = data[3]
        core_inferred_from_sigma[core_inferred_from_sigma=='1'] = '?'
        core_inferred_from_sigma[core_inferred_from_sigma=='0'] = ' '
        distance = data[4].astype(np.float)
        mass_BH = data[5].astype(np.float)
        perr_mass_BH = data[6].astype(np.float)
        merr_mass_BH = data[7].astype(np.float)
	mag_sph = data[8].astype(np.float)
	perr_mag_sph = data[9].astype(np.float)
	merr_mag_sph = data[10].astype(np.float)
	mag_tot = data[11].astype(np.float)
	color = data[12].astype(np.float)
	bar = data[13].astype(np.int)
	delta = data[14].astype(np.float)
		
	## error from 0 to 0.3 mag according to Delta_RMS of profile fit
	## has mean error = 0.08
	goodness = delta - min(delta)
	goodness = goodness/max(goodness)
	err_mag_tot = 0.3*goodness
	       
	log_ML = 3.98*color+0.13 # meidt+2014
	ML = 10**log_ML
	
	mass_sph = ML*10**(-0.4*(mag_sph-3.25))
	perr_mass_sph = (ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) - mass_sph ) 
	merr_mass_sph = (mass_sph - ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) ) 
		
        mmsampletableFile = open(path_tables + 'mmtable_richard.dat', 'w')
        sys.stdout = mmsampletableFile
        
        print '# galaxy  type  core  distance  mass_BH  merr_mass_BH  perr_mass_BH  mag_sph  merr_mag_sph  ',
	print ' perr_mag_sph  mag_tot  err_mag_tot  color  mass_sph  merr_mass_sph  perr_mass_sph '
        for i in range(rows):
        	print gal_name[i], 
        	if bar[i] == 0 or morphtype[i] == 'merger':
        		print morphtype[i], 
        	elif bar[i] == 1:
        		print morphtype[i] + '-bar', 	
        	print str(core[i])+str(core_inferred_from_sigma[i]),  
        	print str("{0:.1f}".format(distance[i])), 
		print mass_BH[i],merr_mass_BH[i],perr_mass_BH[i],
        	print str("{0:.2f}".format(mag_sph[i])), str("{0:.2f}".format(merr_mag_sph[i])), str("{0:.2f}".format(perr_mag_sph[i])),
        	if gal_name[i] in ['m94', 'n3079', 'n4388', 'n4945']:
        		print str("{0:.2f}".format(mag_tot[i])), ' <= ', 
        	else:	
        		print str("{0:.2f}".format(mag_tot[i])), ' 0.25 ',  
        	print str("{0:.2f}".format(color[i])), 
        	print mass_sph[i],merr_mass_sph[i],perr_mass_sph[i]        	
      
        mmsampletableFile.close() 
        terminal = sys.stdout
	
	#print mag_sph[gal_id=='ic1459'],merr_mag_sph[gal_id=='ic1459'],perr_mag_sph[gal_id=='ic1459'],color[gal_id=='ic1459'],mass_sph[gal_id=='ic1459']/10**10,merr_mass_sph[gal_id=='ic1459']/10**10,perr_mass_sph[gal_id=='ic1459']/10**10
        #("{0:.2f}".format(b))
        #("{0:.2f}".format(b))

	


def make_alister_table1():
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, \
		physres.log_r_e_eq_moffat_comb, errV.perr_log_r_e, errV.merr_log_r_e, \
		physres.mu_e_eq_moffat_comb, errV.perr_mu_e, errV.merr_mu_e, \
		physres.log_n_eq_moffat_comb, errV.perr_log_n, errV.merr_log_n, \
		anc.distance \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	log_r_e = data[1].astype(np.float)
	perr_log_r_e = data[2].astype(np.float)
	merr_log_r_e = data[3].astype(np.float)
	mu_e = data[4].astype(np.float)	
	perr_mu_e = data[5].astype(np.float)		
	merr_mu_e = data[6].astype(np.float)		
	log_n = data[7].astype(np.float)			
	perr_log_n = data[8].astype(np.float)	
	merr_log_n = data[9].astype(np.float)	
	distance = data[10].astype(np.float)	
	
        fName = path_tables + 'r_e-mu_e-n.dat'
        f = open(fName, 'w')
        f.write('# galaxy   distance [Mpc]   log_r_e [log(kpc)]   perr_log_r_e   merr_log_r_e   mu_e [mag/arcsec^2]   perr_mu_e   merr_mu_e   log_n   perr_log_n   merr_log_n \n')
        for g, d, r, pr, mr, m, pm, mm, n, pn, mn in zip(gal_id, distance, log_r_e, perr_log_r_e, merr_log_r_e, mu_e, perr_mu_e, merr_mu_e, log_n, perr_log_n, merr_log_n):
        	f.write(str(g) + '     ' + str(d) + '     ' + str(r) + '     ' + str(pr) + '     ' + str(mr) + '     ' + str(m) + '     ' + str(pm) + '     ' + str(mm) + '     ' + str(n) + '     ' + str(pn) + '     ' + str(mn) + '\n')
        f.close()
       
        ############################

def make_shankar_table():
        connection = sql3.connect(dbname)
        cur = connection.cursor()
       
        getdata_query = 'SELECT anc.gal_id, \
		physres.log_r_e_maj_moffat_comb, physres.log_r_e_eq_moffat_comb, errV.perr_log_r_e, errV.merr_log_r_e \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'
	
	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()
	#print data
	gal_id = data[0]
	log_r_e_maj = data[1].astype(np.float)
	log_r_e_eq = data[2].astype(np.float)
	perr_log_r_e = data[3].astype(np.float)
	merr_log_r_e = data[4].astype(np.float)
	
        fName = path_tables + 'r_e_maj-r_e_eq.dat'
        f = open(fName, 'w')
        f.write('# galaxy   log_r_e_maj [log(kpc)]   log_r_e_eq [log(kpc)]   perr_log_r_e   merr_log_r_e   \n')
        for g, a,b,c,d in zip(gal_id, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e):
        	f.write(str(g) + '     ' + str(a) + '     ' + str(b) + '     ' + str(c) + '     ' + str(d)  + '\n')
        f.close()
       
        ############################
	

def make_BHFP_data_table():

	connection = sql3.connect(dbname)
	cur = connection.cursor()

	getdata_query = 'SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.bar, \
		anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, anc.sigma, \
		physres.mag_sph_eq_moffat_comb, physres.mag_tot_eq_moffat_comb, \
		physres.log_n_maj_moffat_comb, physres.log_n_eq_moffat_comb, \
		physres.mu_e_maj_moffat_comb, physres.mu_e_eq_moffat_comb, \
		physres.log_r_e_maj_moffat_comb, physres.log_r_e_eq_moffat_comb, \
		errV.perr_mag_sph, errV.merr_mag_sph, \
		errV.perr_log_n, errV.merr_log_n, \
		errV.perr_mu_e, errV.merr_mu_e, \
	        errV.perr_log_r_e, errV.merr_log_r_e, \
		errV.perr_mu_0, errV.merr_mu_0 \
		FROM Ancillary AS anc \
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id \
		JOIN OneDFitResults AS res ON anc.gal_id = res.gal_id \
		JOIN ErrorsVote as errV ON anc.gal_id = errV.gal_id \
		WHERE anc.fit1D_done = 1;'

	cur.execute(getdata_query)
	datalist = cur.fetchall()
	data= np.asarray(datalist).transpose()

	gal_id = data[0]
	simplemorphtype = data[1]
	core = data[2].astype(np.int)
	bar = data[3].astype(np.int)

	earlytype = np.zeros(len(gal_id))
	for i in range(len(gal_id)):
	    if simplemorphtype[i]=='E' or simplemorphtype[i]=='E/S0' or simplemorphtype[i]=='S0':
	        earlytype[i] = 1
	latetype = np.zeros(len(gal_id))
	for i in range(len(gal_id)):
	    if simplemorphtype[i]=='Sp':
	        latetype[i] = 1

	mbh = data[4].astype(np.float)
	log_mbh = np.log10(mbh)
	perr_mbh = data[5].astype(np.float)
	merr_mbh = data[6].astype(np.float)
	perr_log_mbh = np.log10(1 + perr_mbh/mbh)
	merr_log_mbh = -np.log10(1 - merr_mbh/mbh)

	sigma = data[7].astype(np.float)
	# assign value to n3079
	sigma[gal_id=='n3079'] = 105
	log_sigma = np.log10(sigma)
	err_log_sigma = sigma*[0.0] + np.log10(1.05)

	mag_sph = data[8].astype(np.float)
	mag_tot = data[9].astype(np.float)
	perr_mag_sph = data[16].astype(np.float)
	merr_mag_sph = data[17].astype(np.float)
	err_mag_tot = np.zeros(len(gal_id)) + 0.25

	log_n_maj = data[10].astype(np.float)
	n_maj = 10**log_n_maj
	log_n_eq = data[11].astype(np.float)
	n_eq = 10**log_n_eq
	perr_log_n = data[18].astype(np.float)
	merr_log_n = data[19].astype(np.float)

	mu_e_maj = data[12].astype(np.float)
	mu_e_eq = data[13].astype(np.float)
	perr_mu_e = data[20].astype(np.float)
	merr_mu_e = data[21].astype(np.float)

	# compute mu_0
	b_maj = np.zeros(len(gal_id))
	for i in range(len(b_maj)):
		b_maj[i] = b_n.computeb_n(n_maj[i])
	b_eq = np.zeros(len(gal_id))
	for i in range(len(b_eq)):
		b_eq[i] = b_n.computeb_n(n_eq[i])

	mu_0_maj = mu_e_maj - 2.5*b_maj/np.log(10)
	mu_0_eq = mu_e_eq - 2.5*b_eq/np.log(10)
	perr_mu_0 = data[24].astype(np.float)
	merr_mu_0 = data[25].astype(np.float)

	log_r_e_maj = data[14].astype(np.float)
	#r_e_maj = 10**log_r_e_maj
	log_r_e_eq = data[15].astype(np.float)
	#r_e_eq = 10**log_r_e_eq
	perr_log_r_e = data[22].astype(np.float)
	merr_log_r_e = data[23].astype(np.float)

		
        fiName = '/Users/gsavorgnan/galaxy_vivisection/python_codes/BHFP/BHFP_all_data.dat'
        fi = open(fiName, 'w')
        fi.write('#    log(Mbh)  perr  merr  log(sigma)  err  mag_sph  perr  merr  mag_tot  err  \
		log(n_maj)  log(n_eq)  perr  merr  log(Re_maj)  log(Re_eq)  perr  merr  mu_e_maj  mu_e_eq  perr  merr \
		mu_0_maj  mu_0_eq  perr  merr \n')
        for a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z in zip(log_mbh, perr_log_mbh, merr_log_mbh, 
        	log_sigma, err_log_sigma, mag_sph, perr_mag_sph, merr_mag_sph, mag_tot, err_mag_tot, log_n_maj, log_n_eq, 
        	perr_log_n, merr_log_n, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e, mu_e_maj, mu_e_eq, 
        	perr_mu_e, merr_mu_e, mu_0_maj, mu_0_eq, perr_mu_0, merr_mu_0):
		        	
        	fi.write(str(a) + '     ' + str(b) + '	  ' + str(c) + '     ' + str(d) + '	' + str(e) + '     ' + str(f) )
        	fi.write('     ' + str(g) + '	 ' + str(h) + '     ' + str(i) + '     ' + str(j) + '	  ' + str(k) )
        	fi.write('     ' + str(l) + '	 ' + str(m) + '     ' + str(n) + '     ' + str(o) + '	  ' + str(p) )
        	fi.write('     ' + str(q) + '	 ' + str(r) + '     ' + str(s) + '     ' + str(t) + '	  ' + str(u) )
        	fi.write('     ' + str(v) + '	 ' + str(w) + '     ' + str(x) + '     ' + str(y) + '	  ' + str(z) )
        	fi.write('     ' + '\n')
        fi.close()


        fiName = '/Users/gsavorgnan/galaxy_vivisection/python_codes/BHFP/BHFP_early_data.dat'
        fi = open(fiName, 'w')
        fi.write('#    log(Mbh)  perr  merr  log(sigma)  err  mag_sph  perr  merr  mag_tot  err  \
		log(n_maj)  log(n_eq)  perr  merr  log(Re_maj)  log(Re_eq)  perr  merr  mu_e_maj  mu_e_eq  perr  merr \
		mu_0_maj  mu_0_eq  perr  merr \n')
        for a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,early in zip(log_mbh, perr_log_mbh, merr_log_mbh, 
        	log_sigma, err_log_sigma, mag_sph, perr_mag_sph, merr_mag_sph, mag_tot, err_mag_tot, log_n_maj, log_n_eq, 
        	perr_log_n, merr_log_n, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e, mu_e_maj, mu_e_eq, 
        	perr_mu_e, merr_mu_e, mu_0_maj, mu_0_eq, perr_mu_0, merr_mu_0, earlytype):
		
		if early == 1:
        		fi.write(str(a) + '     ' + str(b) + '	  ' + str(c) + '     ' + str(d) + '	' + str(e) + '     ' + str(f) )
        		fi.write('     ' + str(g) + '	 ' + str(h) + '     ' + str(i) + '     ' + str(j) + '	  ' + str(k) )
        		fi.write('     ' + str(l) + '	 ' + str(m) + '     ' + str(n) + '     ' + str(o) + '	  ' + str(p) )
        		fi.write('     ' + str(q) + '	 ' + str(r) + '     ' + str(s) + '     ' + str(t) + '	  ' + str(u) )
        		fi.write('     ' + str(v) + '	 ' + str(w) + '     ' + str(x) + '     ' + str(y) + '	  ' + str(z) )
        		fi.write('     ' + '\n')
        fi.close()

        fiName = '/Users/gsavorgnan/galaxy_vivisection/python_codes/BHFP/BHFP_ELL_data.dat'
        fi = open(fiName, 'w')
        fi.write('#    log(Mbh)  perr  merr  log(sigma)  err  mag_sph  perr  merr  mag_tot  err  \
		log(n_maj)  log(n_eq)  perr  merr  log(Re_maj)  log(Re_eq)  perr  merr  mu_e_maj  mu_e_eq  perr  merr \
		mu_0_maj  mu_0_eq  perr  merr \n')
        for a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ty in zip(log_mbh, perr_log_mbh, merr_log_mbh, 
        	log_sigma, err_log_sigma, mag_sph, perr_mag_sph, merr_mag_sph, mag_tot, err_mag_tot, log_n_maj, log_n_eq, 
        	perr_log_n, merr_log_n, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e, mu_e_maj, mu_e_eq, 
        	perr_mu_e, merr_mu_e, mu_0_maj, mu_0_eq, perr_mu_0, merr_mu_0, simplemorphtype):
		
		if ty == 'E':
        		fi.write(str(a) + '     ' + str(b) + '	  ' + str(c) + '     ' + str(d) + '	' + str(e) + '     ' + str(f) )
        		fi.write('     ' + str(g) + '	 ' + str(h) + '     ' + str(i) + '     ' + str(j) + '	  ' + str(k) )
        		fi.write('     ' + str(l) + '	 ' + str(m) + '     ' + str(n) + '     ' + str(o) + '	  ' + str(p) )
        		fi.write('     ' + str(q) + '	 ' + str(r) + '     ' + str(s) + '     ' + str(t) + '	  ' + str(u) )
        		fi.write('     ' + str(v) + '	 ' + str(w) + '     ' + str(x) + '     ' + str(y) + '	  ' + str(z) )
        		fi.write('     ' + '\n')
        fi.close()

	
#make_alister_table1()	
#make_BHFP_data_table()	
make_shankar_table()
#make_richard_table()
