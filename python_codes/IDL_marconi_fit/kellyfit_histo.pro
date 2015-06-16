;data = ['mbh_vs_mag_gal_all.dat', 'mbh_vs_mag_gal_early.dat', 'mbh_vs_mag_sph_all.dat', 'mbh_vs_mag_sph_core.dat', 'mbh_vs_mag_sph_sersic.dat', 'mbh_vs_mag_sph_early.dat', 'mbh_vs_mag_sph_late.dat', 'mbh_vs_mag_sph_classbul.dat', 'mbh_vs_mass_sph_early.dat', 'mbh_vs_mass_sph_late.dat']
;data = ['mbh_vs_mag_gal_all.dat', 'mbh_vs_mag_gal_early.dat']
;data = ['mbh_vs_mass_sph_core.dat', 'mbh_vs_mass_sph_sersic.dat']
data = ['mbh_vs_mag_sph_sersic_noout.dat']

openw, 5, 'kellyfit_histo.out'

FOR I = 0,(N_ELEMENTS(data)-1) DO BEGIN & $

	filedata = data[I] & $
	printf, 5, filedata & $
	print, filedata & $
	
	printf, 5, 'y|x' & $
	readcol, filedata, x, xerr, y, yerr & $
	linmix_err, x, y, fit, xsig=xerr, ysig=yerr & $

        cgHistoplot, fit.alpha, BINSIZE=0.01, Output=filedata + '_yx_alpha.ps' & $
        cgHistoplot, fit.beta, BINSIZE=0.01, Output=filedata + '_yx_beta.ps' & $
        cgHistoplot, sqrt(fit.sigsqr), BINSIZE=0.01, Output=filedata + '_yx_epsilon.ps' & $

	a_yx = median(fit.alpha) & $
        err_a_yx = stddev(fit.alpha) & $
        b_yx = median(fit.beta) & $
        err_b_yx = stddev(fit.beta) & $
        e_yx = median(sqrt(fit.sigsqr)) & $
        err_e_yx = stddev(sqrt(fit.sigsqr)) & $

	printf, 5, 'a = ', a_yx, ' \pm', err_a_yx & $
	printf, 5, 'b = ', b_yx, ' \pm', err_b_yx & $
	printf, 5, 'epsilon = ', e_yx, ' \pm', err_e_yx & $
	ss = total((y - a_yx - b_yx*x)^2) & $
	Delta = sqrt(ss/(N_ELEMENTS(x)-2)) & $
	printf, 5, 'Delta = ', Delta & $
	
	printf, 5, '----------------------------------------------' & $

	printf, 5, 'x|y' & $
	linmix_err, (y-mean(y)), x, fit, xsig=yerr, ysig=xerr & $

        cgHistoplot, (fit.alpha), BINSIZE=0.01, Output=filedata + '_xy_alpha.ps' & $
        cgHistoplot, (fit.beta), BINSIZE=0.01, Output=filedata + '_xy_beta.ps' & $
        cgHistoplot, (sqrt(fit.sigsqr)), BINSIZE=0.01, Output=filedata + '_xy_epsilon.ps' & $

	A_inv = median(fit.alpha) & $
        B_inv = median(fit.beta) & $
        E_inv = median(sqrt(fit.sigsqr)) & $
        err_A_inv = stddev(fit.alpha) & $
        err_B_inv = stddev(fit.beta) & $ 
        err_E_inv = stddev(sqrt(fit.sigsqr)) & $
        
	a_xy = -(A_inv/B_inv)+mean(y) & $
        err_a_xy = (err_A_inv^2./B_inv^2. + A_inv^2.*err_B_inv^2./B_inv^4.)^0.5 & $
	b_xy = 1/B_inv & $
	err_b_xy = err_B_inv/B_inv^2 & $
        e_xy = E_inv/(B_inv^2)^0.5 & $  
        err_e_xy = (err_E_inv^2/B_inv^2 + E_inv^2*err_B_inv^2/B_inv^4 )^0.5 & $

	printf, 5, 'a = ', a_xy, ' \pm', err_a_xy & $
	printf, 5, 'b = ', b_xy, ' \pm', err_b_xy & $
	printf, 5, 'epsilon = ', e_xy, ' \pm',  err_e_xy & $
	ss = total((y - a_xy - b_xy*x)^2) & $
	Delta = sqrt(ss/(N_ELEMENTS(x)-2)) & $
	printf, 5, 'Delta = ', Delta & $
	
	printf, 5, '----------------------------------------------' & $

	printf, 5, 'bisector' & $
	b_bis = tan(0.5*(atan(b_yx)+atan(b_xy))) & $
	a_bis = (a_xy-a_yx)*(b_yx-b_bis)/(b_yx-b_xy) + a_yx & $
	err_a_bis = (err_a_yx^2+err_a_xy^2)^0.5/2^0.5 & $
	err_b_bis = tan(err_a_bis+atan(b_bis)) - b_bis & $
	printf, 5, 'a = ', a_bis, ' \pm', err_a_bis & $
	printf, 5, 'b = ', b_bis, ' \pm', err_b_bis & $
	ss = total((y - a_bis - b_bis*x)^2) & $
	Delta = sqrt(ss/(N_ELEMENTS(x)-2)) & $
	printf, 5, 'Delta = ', Delta & $

	printf, 5, '----------------------------------------------' & $


ENDFOR


close,5
