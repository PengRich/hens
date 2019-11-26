if(mod(i, print_number) == 0) then
    k = 0
    do j=1, n_hex
        if(xmin(j, 1)>1.d-3) k = k+1
    enddo
    if(debug==1) then
        write(n_log_file, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        write(n_log_file, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        write(n_log_file, *) i, "p", prob_strategy, sum(prob_strategy)
        write(n_log_file, *) i, "crm", crm
        write(n_log_file, *) i, "crs", crs
        write(n_log_file, *) i, "cfm", cfm
        write(n_log_file, *) i, "cfs", cfs
        write(*, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        write(*, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        write(*, *) i, "p", prob_strategy, sum(prob_strategy)
        write(*, *) i, "crm", crm
        write(*, *) i, "crs", crs
        write(*, *) i, "cfm", cfm
        write(*, *) i, "cfs", cfs
        write(*, *) i, "n_cfr", n_cfr 
    else
        write(n_log_file, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
            sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
        write(*, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
           sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
    endif
endif
