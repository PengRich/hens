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
write(*, '("Time =", f10.1)') finish-start
write(n_log_file, '("ymin =", f15.3)') ymin
write(n_log_file, '("Time =", f10.1)') finish-start
