pro fit_and_print_results, x0, x1, x2, y, x0err, x1err, x2err, yerr
; This function performs Y|(X0,X1) fit and prints the results
    
    x = [[x0-mean(x0)],[x1-mean(x1)],[x2-mean(x2)]]
    
    var_yerr = yerr^2
    
    cov_xerr = MAKE_ARRAY(size(x0, /N_ELEMENTS), 3, 3, /FLOAT)
    
    FOR I = 0, size(x0, /N_ELEMENTS) - 1 DO BEGIN
       sigmaI = [[(x0err[I])^2, 0, 0], [0, (x1err[I])^2, 0], [0, 0, (x2err[I])^2]] & $ 
       cov_xerr[I,*,*] = sigmaI & $ 
    ENDFOR
    
    mlinmix_err, x, y-mean(y), fit, xvar=cov_xerr, yvar=var_yerr


    printf, 6, 'Y|X'
    printf, 6, 'Y - <Y> = a + b0*(X0 - <X0>) + b1*(X1 - <X1>) + b2*(X2 - <X2>) + epsilon'
    a_yx = median(fit.alpha) 
    err_a_yx = stddev(fit.alpha) 
    b0_yx = median(fit.beta[0]) 				   
    err_b0_yx = stddev(fit.beta[0]) 
    b1_yx = median(fit.beta[1]) 				   
    err_b1_yx = stddev(fit.beta[1]) 
    b2_yx = median(fit.beta[2]) 				   
    err_b2_yx = stddev(fit.beta[2]) 
    e_yx = median(sqrt(fit.sigsqr)) 
    err_e_yx = stddev(sqrt(fit.sigsqr)) 
   
    printf, 6, '<X0> = ', mean(x0)
    printf, 6, '<X1> = ', mean(x1)
    printf, 6, '<X2> = ', mean(x2)
    printf, 6, '<Y> = ', mean(y)
    printf, 6, 'a = ', a_yx, ' \pm', err_a_yx 
    printf, 6, 'b0 = ', b0_yx, ' \pm', err_b0_yx 
    printf, 6, 'b1 = ', b1_yx, ' \pm', err_b1_yx 
    printf, 6, 'b2 = ', b2_yx, ' \pm', err_b2_yx 
    printf, 6, 'epsilon = ', e_yx, ' \pm', err_e_yx 
    ss = total(((y-mean(y)) - a_yx - b0_yx*(x0-mean(x0)) - b1_yx*(x1-mean(x1)))^2) 
    Delta = sqrt(ss/(N_ELEMENTS(x1)-4)) 
    printf, 6, 'Delta = ', Delta 
   
    printf, 6, '----------------------------------------------' 
    printf, 6, ' '


end



pro Main3, datafilename, outputfilename

    readcol, datafilename, log_mbh, perr_log_mbh, merr_log_mbh, log_sigma, err_log_sigma, mag_sph, perr_mag_sph, merr_mag_sph, mag_tot, err_mag_tot, log_n_maj, log_n_eq, perr_log_n, merr_log_n, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e, mu_e_maj, mu_e_eq, perr_mu_e, merr_mu_e, mu_0_maj, mu_0_eq, perr_mu_0, merr_mu_0
    
    openw, 6, outputfilename
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_maj - log_n_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, log_n_maj, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_log_n+merr_log_n), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_eq - log_n_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, log_n_eq, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_log_n+merr_log_n), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_maj - mu_e_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, mu_e_maj, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_eq - mu_e_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mu_e_eq, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_maj - mu_0_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, mu_0_maj, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_eq - mu_0_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mu_0_eq, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_maj - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_eq - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_maj - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_maj, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_r_e_eq - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_maj - mu_e_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_maj, mu_e_maj, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_eq - mu_e_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mu_e_eq, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_maj - mu_0_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_maj, mu_0_maj, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_eq - mu_0_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_eq, mu_0_eq, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_maj - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_maj, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_eq - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_eq, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_maj - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_maj, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - log_n_eq - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_n_eq, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_log_n+merr_log_n), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_maj - mu_0_maj'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_maj, mu_0_maj, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_eq - mu_0_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_eq, mu_0_eq, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_maj - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_maj, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_eq - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_eq, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_maj - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_maj, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_e_eq - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_e_eq, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_mu_e+merr_mu_e), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_0_maj - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_0_maj, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_0_eq - mag_sph'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_0_eq, mag_sph, log_mbh, err_log_sigma, 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_0_maj - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_0_maj, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_mu_0+merr_mu_0), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
    printf, 6, 'log_mbh - log_sigma - mu_0_eq - mag_tot'
    printf, 6, ' '
    fit_and_print_results, log_sigma, mu_0_eq, mag_tot, log_mbh, err_log_sigma, 0.5*(perr_mu_0+merr_mu_0), err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)
   
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 6, 'mag_sph - log_sigma - log_r_e_eq - mu_0_eq'
    printf, 6, ' '
    fit_and_print_results, log_sigma, log_r_e_eq, mu_0_eq, mag_sph, err_log_sigma, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_mag_sph+merr_mag_sph)
   
    
    
   
    close, 6

end
