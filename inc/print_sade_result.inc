j = 0
do i=1, n_hex
    if(abs(xmin(i, 1))>1.d-3) then
        j = j + 1
        write(n_log_file, '("Qs[", i3, "]=", f13.3)') i, xmin(i, 1)     
    endif
enddo
n_hex_global = j
y_min_global = ymin
write(*, '("ymin =", f15.3)') ymin
write(n_log_file, '("ymin =", f15.3)') ymin
write(*, *) "p", prob_strategy, sum(prob_strategy)
write(*, *) "crm", crm
write(*, *) "crs", crs
write(*, *) "cfm", cfm
write(*, *) "cfs", cfs
