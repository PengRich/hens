if(mod(i, print_number) == 0) then
    k = 0
    do j=1, n_hex
        if(xmin(j, 1)>1.d-3) k = k+1
    enddo
    if(debug==1) then
        write(13, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        write(13, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        write(13, *) i, "p", prob_strategy, sum(prob_strategy)
        write(13, *) i, "crm", crm
        write(13, *) i, "crs", crs
        write(13, *) i, "cfm", cfm
        write(13, *) i, "cfs", cfs
        write(*, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        write(*, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        write(*, *) i, "p", prob_strategy, sum(prob_strategy)
        write(*, *) i, "crm", crm
        write(*, *) i, "crs", crs
        write(*, *) i, "cfm", cfm
        write(*, *) i, "cfs", cfs
        write(*, *) i, "n_cfr", n_cfr 
    else
        write(13, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
            sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
        write(*, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
           sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
    endif
endif
