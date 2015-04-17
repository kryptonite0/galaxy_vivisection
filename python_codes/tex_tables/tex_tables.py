import sqlite3 as sql3
import os
import sys
import numpy as np

terminal = sys.stdout
dbname = '/Users/gsavorgnan/galaxy_vivisection/python_codes/databases/galaxy_vivisection.db'
sampletableFileName = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/table_sample.tex'
fitresultstableFileName = '/Users/gsavorgnan/galaxy_vivisection/papers/data_paper/table_fitresults.tex'
mmsampletableFileName = '/Users/gsavorgnan/galaxy_vivisection/papers/MbhMsph/table_sample.tex'

def putBlankInPlaceOfNone(entry):
        if entry is None:
                entry = ' '
        return entry  
	
def header_datapaper_sampletable(caption):

        print '\\begin{table*}                                        '
        print '\\small                                                '
        print '\\begin{center}                                        '
	if caption:
        	print '\\caption{{\\bf Galaxy sample.}                        '
        	print '\\emph{Column (1):} Galaxy name.                       '
        	print '\\emph{Column (2):} Distance.                                   '
        	print '\\emph{Column (3):} Black hole mass.                                   '
       		print '\\emph{Column (4):} Reference of the black hole mass reported here \
			(G+03 = \citealt{greenhill2003}, GS13 = \citealt{grahamscott2013}; R+13b = \citealt{rusli2013bhmassesDM}).                                   '
        	print '\\emph{Column (5):} Presence of a partially depleted core. \
			The question mark is used when the classification has come from the velocity dispersion criteria mentioned in Section \\ref{sec:corser}. \
			The value of the core break radius is reported in parenthesis when available.  '
        	print '\\emph{Column (6):} Reference of the identification of a partially depleted core \
			(G+94 = \citealt{grillmair1994}; \
			F+97 = \citealt{forbes1997}; \
			Q+00 = \citealt{quillen2000}, \
			T+04 = \citealt{trujillo2004coresersicmodel}; \
			F+06 = \citealt{ferrarese2006acsvcs}; \
			J+11 = \citealt{jardel2011}; \
			R+11 = \citealt{richings2011}; \
			R+13a = \citealt{rusli2013}).  '
		print '\\emph{Column (7):} Kinematical classification (fast/slow rotator).'
		print '\\emph{Column (8):} Availability of velocity map (A = ATLAS$^{\\rm 3D}$, S = SLUGGS). '
		print '\\emph{Column (9):} Completion of 1D fit. '
		print '\\emph{Column (10):} Completion of 2D fit. }                                 '
        print '\\begin{tabular}{llllllllll}                           '
        print '\\hline                                                '
        print '\\multicolumn{1}{l}{{\\bf Galaxy}} &                   '
        print '\\multicolumn{1}{l}{{\\bf Distance}} &                 '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{M_{\\rm BH}}$}} &  '
        print '\\multicolumn{1}{l}{{\\bf Ref.}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Core}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Ref.}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Rot.}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Vel. map}} &                 '
        print '\\multicolumn{1}{l}{{\\bf 1D fit}} &                   '
        print '\\multicolumn{1}{l}{{\\bf 2D fit}} \\\\                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{[Mpc]} &                           '
        print '\\multicolumn{1}{l}{$[10^8~\\rm M_{\odot}]$} &         '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} \\\\                             '
        print '\\multicolumn{1}{l}{(1)} &                             '
        print '\\multicolumn{1}{l}{(2)} &                             '
        print '\\multicolumn{1}{l}{(3)} &                             '
        print '\\multicolumn{1}{l}{(4)} &                             '
        print '\\multicolumn{1}{l}{(5)} &                             '
        print '\\multicolumn{1}{l}{(6)} &                             '
        print '\\multicolumn{1}{l}{(7)} &                             '
        print '\\multicolumn{1}{l}{(8)} &                             '
        print '\\multicolumn{1}{l}{(9)} &                             '
        print '\\multicolumn{1}{l}{(10)} \\\\                         '
        print '\\hline                                                '

def footer_datapaper_sampletable(label):

        print '\\hline         '
        print '\\end{tabular}   '
        if label:
		print '\\label{tab:sample} '
        print '\\end{center}    '
        print '\\end{table*}    '

def header_mmpaper_sampletable(caption):

        print '\\begin{table*}                                        '
        print '\\small                                                '
        print '\\begin{center}                                        '
	if caption:
        	print '\\caption{{\\bf Galaxy sample.}                        '
        	print '\\emph{Column (1):} Galaxy name.                       '
		print '\\emph{Column (2):} Morphological type (E=elliptical, S0(B)=(barred) lenticular, Sp(B)=(barred) spiral, merger=merger).                       '
        	print '\\emph{Column (3):} Presence of a partially depleted core. \
			The question mark is used when the classification has come from the velocity dispersion criteria mentioned in Section \\ref{sec:data}.   '
        	print '\\emph{Column (4):} Distance.                                   '
        	print '\\emph{Column (5):} Black hole mass.                                   '
		print '\\emph{Column (6):} Absolute $3.6\\rm~\mu m$ bulge magnitude.                                   '
		print '\\emph{Column (7):} Absolute $3.6\\rm~\mu m$ galaxy magnitude. \
			The four galaxy magnitudes marked with a * are upper limits.                                   '
		print '\\emph{Column (8):} $[3.6]-[4.5]$ colour.                                   '
		print '\\emph{Column (9):} Bulge stellar mass. }                      '
        print '\\begin{tabular}{lllllllll}                           '
        print '\\hline                                                '
        print '\\multicolumn{1}{l}{{\\bf Galaxy}} &                   '
        print '\\multicolumn{1}{l}{{\\bf Type}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Core}} &                     '
        print '\\multicolumn{1}{l}{{\\bf Distance}} &                 '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{M_{\\rm BH}}$}} &  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{MAG_{\\rm sph}}$}} &  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{MAG_{\\rm gal}}$}} &  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{[3.6]-[4.5]}$}} &  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{M_{\\rm *,sph}}$}} \\\\  '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{} &                                '
        print '\\multicolumn{1}{l}{[Mpc]} &                           '
        print '\\multicolumn{1}{l}{$[10^8~\\rm M_{\odot}]$} &         '
        print '\\multicolumn{1}{l}{[mag]} &                                '
        print '\\multicolumn{1}{l}{[mag]} &                                '
        print '\\multicolumn{1}{l}{[mag]} &                                '
        print '\\multicolumn{1}{l}{$[10^{10}~\\rm M_{\odot}]$} \\\\                             '
        print '\\multicolumn{1}{l}{(1)} &                             '
        print '\\multicolumn{1}{l}{(2)} &                             '
        print '\\multicolumn{1}{l}{(3)} &                             '
        print '\\multicolumn{1}{l}{(4)} &                             '
        print '\\multicolumn{1}{l}{(5)} &                             '
        print '\\multicolumn{1}{l}{(6)} &                             '
        print '\\multicolumn{1}{l}{(7)} &                             '
        print '\\multicolumn{1}{l}{(8)} &                             '
        print '\\multicolumn{1}{l}{(9)} \\\\                         '
        print '\\hline                                                '

def footer_mmpaper_sampletable(label):

        print '\\hline         '
        print '\\end{tabular}   '
        if label:
		print '\\label{tab:sample} '
        print '\\end{center}    '
        print '\\end{table*}    '

	  

def datapaper_sampletable():
        connection = sql3.connect(dbname)
        cur = connection.cursor()

        cur.execute('''SELECT A.gal_id, A.distance, A.mass_BH, A.perr_mass_BH, A.merr_mass_BH, A.source, 
                A.core, A.core_inferred_from_sigma, A.size_core, A.ref_size_core, A.fit1D_done, A.fit2D_done,
                K.fastrotator, K.slowrotator, K.velmap_atlas3d, K.velmap_slugg
                FROM Ancillary as A 
                JOIN Kinematics as K ON A.gal_id = K.gal_id
                WHERE A.IRAC_image = 1
                ORDER BY A.gal_id''')
                
        datat = cur.fetchall()
        data= np.asarray(datat).transpose()
        
        shape = np.shape(data)
        columns = shape[0]
        rows = shape[1]
        for i in range(columns):
                for j in range(rows):
                        data[i,j] = putBlankInPlaceOfNone(data[i,j])
        
        gal_name = data[0]
        distance = data[1]
        mass_BH = data[2]/10**8
        perr_mass_BH = data[3]/10**8
        merr_mass_BH = data[4]/10**8
        source = data[5]
        core = data[6]
        core[core==1] = 'yes'
        core[core==0] = 'no'
        core_inferred_from_sigma = data[7]
        core_inferred_from_sigma[core_inferred_from_sigma==1] = '?'
        core_inferred_from_sigma[core_inferred_from_sigma==0] = ' '
        size_core = data[8]
        ref_size_core = data[9]
        fit1D_done = data[10]
        fit1D_done[fit1D_done==1] = 'yes'
        fit1D_done[fit1D_done==0] = 'no'
        fit2D_done = data[11]
        fit2D_done[fit2D_done==1] = 'yes'
        fit2D_done[fit2D_done==0] = 'no'
        fastrotator = data[12]
        fastrotator[fastrotator==1] = 'FAST'
        fastrotator[fastrotator==0] = ' '
        slowrotator = data[13]
        slowrotator[slowrotator==1] = 'SLOW'
        slowrotator[slowrotator==0] = ' '
        velmap_atlas3d = data[14]
        velmap_slugg = data[15]
        
        sampletableFile = open(sampletableFileName, 'w')
        sys.stdout = sampletableFile
        
        header_datapaper_sampletable(True)
        for i in range(rows):
		if gal_name[i] == 'n3384':
			footer_datapaper_sampletable(True)
			print
			header_datapaper_sampletable(False)
                gal_name[i] = gal_name[i].replace('circinus', 'Circinus ')
                if gal_name[i][0] == 'n':
                        gal_name[i] = gal_name[i].replace('n', 'NGC ')
                if gal_name[i][0] == 'i':
                        gal_name[i] = gal_name[i].replace('ic', 'IC ')
                if gal_name[i][0] == 'u':
                        gal_name[i] = gal_name[i].replace('ugc', 'UGC ')
                if gal_name[i][0] == 'm':
                        gal_name[i] = gal_name[i].replace('m', 'M')
                if gal_name[i][-1] == 'a':
                        gal_name[i] = gal_name[i].replace('a', 'A')
                gal_name[i] = gal_name[i].replace('exp', '')    
                print gal_name[i], ' & ', 
                print '$'+str("{0:.1f}".format(distance[i]))+'$', ' & ', 
                print_bhmass(mass_BH[i],merr_mass_BH[i],perr_mass_BH[i])
                #print '$'+str("{0:.4f}".format(mass_BH[i]))+'_{-'+str("{0:.4f}".format(merr_mass_BH[i]))+'}^{+'+str("{0:.4f}".format(perr_mass_BH[i]))+'}$', ' & ', 
                print str(source[i]), ' & ', 
                print str(core[i])+str(core_inferred_from_sigma[i]), 
		if str(size_core[i]) is not ' ':
			print '$('+(str(size_core[i])).replace('.','".')+')$', ' & ', 
		else:
			print ' & ',	
		print str(ref_size_core[i].replace('&','')), ' & ',
                print str(fastrotator[i])+str(slowrotator[i]), ' & ',
                if velmap_atlas3d[i] == 1 and velmap_slugg[i] == 1:
                        print 'A, S', ' & ', 
                else:
                        if velmap_atlas3d[i] == 1:
                                print 'A', ' & ',
                        elif velmap_slugg[i] == 1:
                                print 'S', ' & ',
                        else: 
                                print ' ', ' & ',                       
                print fit1D_done[i], ' & ', fit2D_done[i], 
                print ' \\\\ '
        footer_datapaper_sampletable(False)
	
        sampletableFile.close() 
        terminal = sys.stdout
        #("{0:.2f}".format(b))
        #("{0:.2f}".format(b))
        
def print_bhmass(mass,merr,perr):
        mass = str("{0:.5f}".format(mass))
        merr = str("{0:.5f}".format(merr))
        perr = str("{0:.5f}".format(perr))
        for i in range(1,len(merr)):
                if merr[-i] is not '0':
                        break   
        break_index_merr = i
        for i in range(1,len(perr)):
                if perr[-i] is not '0':
                        break   
        break_index_perr = i    
        break_index = min(break_index_merr, break_index_perr)   
        
        mass_trunc = mass[:-break_index+1]
        merr_trunc = merr[:-break_index+1]
        perr_trunc = perr[:-break_index+1]
        
        if mass_trunc[-1] == '.':
                mass_trunc = mass_trunc[:-1]
        if merr_trunc[-1] == '.':
                merr_trunc = merr_trunc[:-1]
        if perr_trunc[-1] == '.':
                perr_trunc = perr_trunc[:-1]
        
        print '$' + mass_trunc + '_{-' + merr_trunc + '}^{+' + perr_trunc + '}$ ', ' & ',

def print_mass(mass,merr,perr):

	if mass<1:
		print '$' + str("{0:.2f}".format(mass)) + '_{' + str("{0:.2f}".format(merr)) + '}^{+' + str("{0:.2f}".format(perr)) + '}$ ', 
	if mass>1 and mass<10:
		print '$' + str("{0:.1f}".format(mass)) + '_{' + str("{0:.1f}".format(merr)) + '}^{+' + str("{0:.1f}".format(perr)) + '}$ ', 
	if mass>10:
		print '$' + str("{0:.0f}".format(mass)) + '_{' + str("{0:.0f}".format(merr)) + '}^{+' + str("{0:.0f}".format(perr)) + '}$ ', 
		
def header_datapaper_fitresultstable(caption):

        print '\\begin{table*}                                        '
        print '\\small                                                '
        print '\\begin{center}                                        '
	if caption:
        	print '\\caption{{\\bf Results of galaxy decompositions.}                        '
        	print '\\emph{Column (1):} Galaxy name.                       '
        	print '\\emph{Column (2-4):} Effective radius (in units of $[\\rm arcsec]$), ',
		print 'surface brightness at the effective radius (in units of $[\\rm mag~arcsec^{-2}]$) ',
		print 'and S\\\'ersic index for 1D fits along the major-axis.'
        	print '\\emph{Column (5-9):} Effective radius (in units of $[\\rm arcsec]$), ',
		print 'surface brightness at the effective radius (in units of $[\\rm mag~arcsec^{-2}]$), ',
		print 'S\\\'ersic index, ', 
		print 'spheroid apparent magnitude (in units of $[\\rm mag]$) ',
		print 'and galaxy apparent magnitude (in units of $[\\rm mag]$) for 1D fits along the equivalent-axis.'
		print '\\emph{Column (10):} Quality flag of the 1D fits (see Section \\ref{sec:err}). '
		print '\\emph{Column (11-13):}  Effective radius (in units of $[\\rm arcsec]$), ',
		print 'S\\\'ersic index, ', 
		print 'and spheroid apparent magnitude (in units of $[\\rm mag]$) for 2D fits. } '
        print '\\begin{tabular}{lllllllllllll}                           '
        print '\\hline                                                '
	print ' & '
	print '\\multicolumn{3}{l}{{\\bf 1D Major-axis}} &                   '
	print '\\multicolumn{5}{l}{{\\bf 1D Equivalent-axis }} &                   '
	print ' & '
	print '\\multicolumn{3}{l}{{\\bf 2D }} \\\\                    '
        print '\\multicolumn{1}{l}{{\\bf Galaxy }} &                   '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{R_{\\rm e}}$ }} &                 '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{\\mu_{\\rm e}}$ }} &  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{n}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{R_{\\rm e}}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{\\mu_{\\rm e}}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{n}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{m_{\\rm sph}}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{m_{\\rm gal}}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf Q.F. }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{R_{\\rm e}}$ }} &			  '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{n}$ }} &                   '
        print '\\multicolumn{1}{l}{{\\bf $\\bm{m_{\\rm sph}}$ }} \\\\                '
       #print '\\multicolumn{1}{l}{} &  			      '
       #print '\\multicolumn{1}{l}{$[\\rm arcsec]$} &				'
       #print '\\multicolumn{1}{l}{$[\\rm mag~arcsec^{-2}]$} &         '
       #print '\\multicolumn{1}{l}{} &  			      '
       #print '\\multicolumn{1}{l}{$[\\rm arcsec]$} &				     '
       #print '\\multicolumn{1}{l}{$[\\rm mag~arcsec^{-2}]$} &  			      '
       #print '\\multicolumn{1}{l}{} &  			      '
       #print '\\multicolumn{1}{l}{$[\\rm mag]$} &				  '
       #print '\\multicolumn{1}{l}{$[\\rm mag]$} &				  '
       #print '\\multicolumn{1}{l}{$[\\rm arcsec]$} &				     '
       #print '\\multicolumn{1}{l}{} &  			      '
       #print '\\multicolumn{1}{l}{$[\\rm mag]$} \\\\				  '
        print '\\multicolumn{1}{l}{(1)} &                             '
        print '\\multicolumn{1}{l}{(2)} &                             '
        print '\\multicolumn{1}{l}{(3)} &                             '
        print '\\multicolumn{1}{l}{(4)} &                             '
        print '\\multicolumn{1}{l}{(5)} &                             '
        print '\\multicolumn{1}{l}{(6)} &                             '
        print '\\multicolumn{1}{l}{(7)} &                             '
        print '\\multicolumn{1}{l}{(8)} &                             '
        print '\\multicolumn{1}{l}{(9)} &                             '
        print '\\multicolumn{1}{l}{(10)} &                             '
        print '\\multicolumn{1}{l}{(11)} &                             '
        print '\\multicolumn{1}{l}{(12)} &                             '
        print '\\multicolumn{1}{l}{(13)} \\\\                         '
        print '\\hline                                                '

def footer_datapaper_fitresultstable(label):

        print '\\hline         '
        print '\\end{tabular}   '
        if label:
		print '\\label{tab:fitres} '
        print '\\end{center}    '
        print '\\end{table*}    '


def datapaper_fitresultstable():
	connection = sql3.connect(dbname)
        cur = connection.cursor()

        cur.execute('''SELECT anc.gal_id, 
		one.r_e_maj_moffat_comb, one.mu_e_maj_moffat_comb, one.n_maj_moffat_comb,
	 	one.r_e_eq_moffat_comb, one.mu_e_eq_moffat_comb, one.n_eq_moffat_comb, one.mag_sph_eq_moffat_comb, one.mag_tot_eq_moffat_comb,
		two.r_e_moffat, two.n_moffat, two.mag_sph_moffat,
		anc.fit1D_done, anc.fit2D_done,
		err.err_n_vote 	
                FROM Ancillary as anc 
                JOIN OneDFitResults as one ON anc.gal_id = one.gal_id
		JOIN TwoDFitResults as two ON anc.gal_id = two.gal_id
		JOIN ErrorsVote as err ON anc.gal_id = err.gal_id
                WHERE anc.IRAC_image = 1
                ORDER BY anc.gal_id''')
                
        datat = cur.fetchall()
        data= np.asarray(datat).transpose()
	
        shape = np.shape(data)
        columns = shape[0]
        rows = shape[1]
        #for i in range(columns):
        #        for j in range(rows):
        #                data[i,j] = putBlankInPlaceOfNone(data[i,j])
        
        gal_name = data[0]
        r_e_maj_1D = data[1]
	mu_e_maj_1D = data[2]
	n_maj_1D = data[3]
	r_e_eq_1D = data[4]
	mu_e_eq_1D = data[5]
	n_eq_1D = data[6]
	mag_sph_eq_1D = data[7]
	mag_tot_eq_1D = data[8]
	r_e_2D = data[9]
	n_2D = data[10]
	mag_sph_2D = data[11]
	fit1D_done = data[12]
	fit2D_done = data[13]
	err_vote = data[14]
	
	fitresultstableFile = open(fitresultstableFileName, 'w')
        sys.stdout = fitresultstableFile
        
        header_datapaper_fitresultstable(True)
        for i in range(rows):
		if gal_name[i] == 'n3384':
			footer_datapaper_fitresultstable(True)
			print
			header_datapaper_fitresultstable(False)
                gal_name[i] = gal_name[i].replace('circinus', 'Circinus ')
                if gal_name[i][0] == 'n':
                        gal_name[i] = gal_name[i].replace('n', 'NGC ')
                if gal_name[i][0] == 'i':
                        gal_name[i] = gal_name[i].replace('ic', 'IC ')
                if gal_name[i][0] == 'u':
                        gal_name[i] = gal_name[i].replace('ugc', 'UGC ')
                if gal_name[i][0] == 'm':
                        gal_name[i] = gal_name[i].replace('m', 'M')
                if gal_name[i][-1] == 'a':
                        gal_name[i] = gal_name[i].replace('a', 'A')
                gal_name[i] = gal_name[i].replace('exp', '')    
		print gal_name[i], ' \\quad & ', 
                if r_e_maj_1D[i] is not None and fit1D_done[i]==1:
                	print '$'+str("{0:.1f}".format(r_e_maj_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if mu_e_maj_1D[i] is not None and fit1D_done[i]==1:
			print '$'+str("{0:.2f}".format(mu_e_maj_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if n_maj_1D[i] is not None and fit1D_done[i]==1:
			print '$'+str("{0:.1f}".format(n_maj_1D[i]))+'$', ' \\quad \\quad & ', 
		else:
			print ' -- \\quad \\quad & ',
		if r_e_eq_1D[i] is not None and fit1D_done[i]==1:	
			print '$'+str("{0:.1f}".format(r_e_eq_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if mu_e_eq_1D[i] is not None and fit1D_done[i]==1:	
			print '$'+str("{0:.2f}".format(mu_e_eq_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if n_eq_1D[i] is not None and fit1D_done[i]==1:		
			print '$'+str("{0:.1f}".format(n_eq_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if mag_sph_eq_1D[i] is not None and fit1D_done[i]==1:
			print '$'+str("{0:.2f}".format(mag_sph_eq_1D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if mag_tot_eq_1D[i] is not None and fit1D_done[i]==1:	
			print '$'+str("{0:.2f}".format(mag_tot_eq_1D[i]))+'$', ' \\quad \\quad & ', 
		else:
			print ' -- \\quad \\quad & '
		if err_vote[i] is not None and fit1D_done[i]==1:
			print '$'+str(err_vote[i])+'$', ' \\quad \\quad & ',	
		else:
			print ' -- \\quad \\quad & '	
		if r_e_2D[i] is not None and fit2D_done[i]==1:
			print '$'+str("{0:.1f}".format(r_e_2D[i]))+'$', ' & ',
		else:
			print ' -- & ',
		if n_2D[i] is not None and fit2D_done[i]==1:
			print '$'+str("{0:.1f}".format(n_2D[i]))+'$', ' & ', 
		else:
			print ' -- & ',
		if mag_sph_2D[i] is not None and fit2D_done[i]==1:
			print '$'+str("{0:.2f}".format(mag_sph_2D[i]))+'$', 
                else:
			print ' --  ',
		print ' \\\\ '
        footer_datapaper_fitresultstable(False)
	
        fitresultstableFile.close() 
        terminal = sys.stdout
        #("{0:.2f}".format(b))
        #("{0:.2f}".format(b))

def mmpaper_sampletable():
	connection = sql3.connect(dbname)
        cur = connection.cursor()

        cur.execute('''SELECT anc.gal_id, anc.simplemorphtype, anc.core, anc.core_inferred_from_sigma, anc.distance,  
		anc.mass_BH, anc.perr_mass_BH, anc.merr_mass_BH, 
		physres.mag_sph_eq_moffat_comb, 
		errV.perr_mag_sph, errV.merr_mag_sph, 
		physres.mag_tot_eq_moffat_comb, 
		col.color, anc.bar
		FROM Ancillary AS anc 
		JOIN OneDFitResultsPhysicalUnits AS physres ON anc.gal_id = physres.gal_id 
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
	morphtype = data[1]
        core = data[2]
        core[core=='1'] = 'yes'
        core[core=='0'] = 'no'
        core_inferred_from_sigma = data[3]
        core_inferred_from_sigma[core_inferred_from_sigma=='1'] = '?'
        core_inferred_from_sigma[core_inferred_from_sigma=='0'] = ' '
        distance = data[4].astype(np.float)
        mass_BH = data[5].astype(np.float)/10**8
        perr_mass_BH = data[6].astype(np.float)/10**8
        merr_mass_BH = data[7].astype(np.float)/10**8
	mag_sph = data[8].astype(np.float)
	perr_mag_sph = data[9].astype(np.float)
	merr_mag_sph = data[10].astype(np.float)
	mag_tot = data[11].astype(np.float)
	color = data[12].astype(np.float)
	bar = data[13].astype(np.int)
		
	log_ML = 3.98*color+0.13 # meidt+2014
	ML = 10**log_ML
	
	mass_sph = ML*10**(-0.4*(mag_sph-3.25))/10**10
	perr_mass_sph = (ML*10**(-0.4*(mag_sph-merr_mag_sph-3.25)) - mass_sph ) /10**10
	merr_mass_sph = (mass_sph - ML*10**(-0.4*(mag_sph+perr_mag_sph-3.25)) ) /10**10
		
	mmsampletableFile = open(mmsampletableFileName, 'w')
        sys.stdout = mmsampletableFile
        
        header_mmpaper_sampletable(True)
        for i in range(rows):
		if gal_name[i] == 'n4473':
			footer_mmpaper_sampletable(True)
			print
			header_mmpaper_sampletable(False)
                gal_name[i] = gal_name[i].replace('circinus', 'Circinus ')
                if gal_name[i][0] == 'n':
                        gal_name[i] = gal_name[i].replace('n', 'NGC ')
                if gal_name[i][0] == 'i':
                        gal_name[i] = gal_name[i].replace('ic', 'IC ')
                if gal_name[i][0] == 'u':
                        gal_name[i] = gal_name[i].replace('ugc', 'UGC ')
                if gal_name[i][0] == 'm':
                        gal_name[i] = gal_name[i].replace('m', 'M')
                if gal_name[i][-1] == 'a':
                        gal_name[i] = gal_name[i].replace('a', 'A')
                gal_name[i] = gal_name[i].replace('exp', '')    
                print gal_name[i], ' & ', 
		if bar[i] == 0 or morphtype[i] == 'merger':
			print morphtype[i], ' & ', 
		elif bar[i] == 1:
			print morphtype[i] + 'B', ' & ', 	
                print str(core[i])+str(core_inferred_from_sigma[i]), ' & ',  
                print '$'+str("{0:.1f}".format(distance[i]))+'$', ' & ', 
                print_bhmass(mass_BH[i],merr_mass_BH[i],perr_mass_BH[i])
		print '$' + str("{0:.2f}".format(mag_sph[i])) + '_{-' + str("{0:.2f}".format(merr_mag_sph[i])) + '}^{+' + str("{0:.2f}".format(perr_mag_sph[i])) + '}$ ', ' & ',
		print '$' + str("{0:.2f}".format(mag_tot[i])) + '$ ', 
		if gal_name[i] in ['M94', 'NGC 3079', 'NGC 4388', 'NGC 4945']:
			print '*',
		print ' & ',
		print '$' + str("{0:.2f}".format(color[i])) +'$', ' & ', 
		#print '$' + str("{0:.2f}".format(mass_sph[i])) + '_{' + str("{0:.2f}".format(merr_mass_sph[i])) + '}^{+' + str("{0:.2f}".format(perr_mass_sph[i])) + '}$ ', 
                print_mass(mass_sph[i],merr_mass_sph[i],perr_mass_sph[i])
		print ' \\\\ '
        footer_mmpaper_sampletable(False)
	
        mmsampletableFile.close() 
        terminal = sys.stdout
        #("{0:.2f}".format(b))
        #("{0:.2f}".format(b))



        
def main():
        #datapaper_sampletable()
	#datapaper_fitresultstable()
	mmpaper_sampletable()

main()
