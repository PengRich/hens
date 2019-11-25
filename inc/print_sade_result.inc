j = 0
do i=1, n_hex
    if(abs(xmin(i, 1))>1.d-3) then
        j = j + 1
        write(13, '("Qs[", i3, "]=", f13.3)') i, xmin(i, 1)     
    endif
enddo
n_hex_global = j
y_min_global = ymin
write(*, '("ymin =", f15.3)') ymin
write(*, '("Time =", f10.1)') finish-start
write(13, '("ymin =", f15.3)') ymin
write(13, '("Time =", f10.1)') finish-start
