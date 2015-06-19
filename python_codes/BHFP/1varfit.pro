pro fit_and_print_results, x, y, xerr, yerr
; This function performs Y|X fit and prints the results

    linmix_err, x-mean(x), y-mean(y), fit, xsig=xerr, ysig=yerr
    printf, 5, 'Y|X'
    printf, 5, 'Y - <Y> = a + b*(X - <X>) + epsilon'
    a_yx = median(fit.alpha) 
    err_a_yx = stddev(fit.alpha) 
    b_yx = median(fit.beta)					   
    err_b_yx = stddev(fit.beta) 
    e_yx = median(sqrt(fit.sigsqr)) 
    err_e_yx = stddev(sqrt(fit.sigsqr)) 

    printf, 5, '<X> = ', mean(x)
    printf, 5, '<Y> = ', mean(y)
    printf, 5, 'a = ', a_yx, ' \pm', err_a_yx 
    printf, 5, 'b = ', b_yx, ' \pm', err_b_yx 
    printf, 5, 'epsilon = ', e_yx, ' \pm', err_e_yx 
    ss = total(((y-mean(y)) - a_yx - b_yx*(x-mean(x)))^2) 
    Delta = sqrt(ss/(N_ELEMENTS(x)-2)) 
    printf, 5, 'Delta = ', Delta 

    printf, 5, '----------------------------------------------' 
    printf, 5, ' '

end

pro Main1, datafilename, outputfilename

    readcol, datafilename, log_mbh, perr_log_mbh, merr_log_mbh, log_sigma, err_log_sigma, mag_sph, perr_mag_sph, merr_mag_sph, mag_tot, err_mag_tot, log_n_maj, log_n_eq, perr_log_n, merr_log_n, log_r_e_maj, log_r_e_eq, perr_log_r_e, merr_log_r_e, mu_e_maj, mu_e_eq, perr_mu_e, merr_mu_e, mu_0_maj, mu_0_eq, perr_mu_0, merr_mu_0

    openw, 5, outputfilename

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - log_sigma'
    printf, 5, ' '
    fit_and_print_results, log_sigma, log_mbh, err_log_sigma, 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mag_sph'
    printf, 5, ' '
    fit_and_print_results, mag_sph, log_mbh, 0.5*(perr_mag_sph+merr_mag_sph), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mag_tot'
    printf, 5, ' '
    fit_and_print_results, mag_tot, log_mbh, err_mag_tot, 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - log_r_e_maj'
    printf, 5, ' '
    fit_and_print_results, log_r_e_maj, log_mbh, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - log_r_e_eq'
    printf, 5, ' '
    fit_and_print_results, log_r_e_eq, log_mbh, 0.5*(perr_log_r_e+merr_log_r_e), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - log_n_maj'
    printf, 5, ' '
    fit_and_print_results, log_n_maj, log_mbh, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - log_n_eq'
    printf, 5, ' '
    fit_and_print_results, log_n_eq, log_mbh, 0.5*(perr_log_n+merr_log_n), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mu_e_maj'
    printf, 5, ' '
    fit_and_print_results, mu_e_maj, log_mbh, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mu_e_eq'
    printf, 5, ' '
    fit_and_print_results, mu_e_eq, log_mbh, 0.5*(perr_mu_e+merr_mu_e), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mu_0_maj'
    printf, 5, ' '
    fit_and_print_results, mu_0_maj, log_mbh, 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    printf, 5, 'log_mbh - mu_0_eq'
    printf, 5, ' '
    fit_and_print_results, mu_0_eq, log_mbh, 0.5*(perr_mu_0+merr_mu_0), 0.5*(perr_log_mbh+merr_log_mbh)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
    close, 5

end
