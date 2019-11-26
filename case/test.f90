program test 
use random_numbers
use sade
implicit none
integer(kind=4), parameter :: n=8
integer(kind=4) :: i, j, k, np, idx_topo(n)
real(kind=8) :: xt, yt 
real(kind=8), allocatable :: Qs(:)

    ! call execute_command_line ("ls")
    print *, "set seed to 1"
    call set_random_seed(1)
    do i=1, 3
        call random_number(xt)
        print *, xt
    enddo
    call set_random_seed(1)
    do i=1, 3
        call random_number(xt)
        print *, xt
    enddo

    print *, "not set seed"
    call set_random_seed()
    do i=1, 3
        call random_number(xt)
        print *, xt
    enddo
    
    ! test get_log_filename
    call get_log_filename("hello")
    print *, log_filename
    
    print *, "start point is 0.1"
    xt = 0.1d0
    do i=1, 5 
        print *, rand(xt)
    enddo
    print *, "start point is 0.1"
    xt = 0.1d0
    do i=1, 5 
        print *, rand(xt)
    enddo
    print *, "start point is 0.2"
    xt = 0.2d0
    do i=1, 5 
        print *, rand(xt)
    enddo
    
    j = 0
    k = 0
    do i=1, 5000
        j = j +1
        if(rand(xt) < 0.3d0) then
            k = k + 1
        endif
    enddo
    print *, "simulated 0.3 is ", real(k)/real(j)

    call init_case("10sp1", 4)
    allocate(Qs(n_hex))
    include "result_10sp1_43392.inc"
    print *, tac(Qs)
    deallocate(Qs)
    call deallocate_var()

    call init_case("9sp1", 4)
    allocate(Qs(n_hex))
    include "result_9sp1_2935m.inc"
    print *, tac(Qs)
    deallocate(Qs)
    call deallocate_var()

    call init_case("9sp1", 4)
    allocate(Qs(n_hex))
    include "result_9sp1_2942m.inc"
    print *, tac(Qs)
    deallocate(Qs)
    call deallocate_var()

    call init_case("10sp2", 4)
    allocate(Qs(n_hex))
    include "result_10sp2_5640616.inc"
    print *, tac(Qs)
    deallocate(Qs)
    call deallocate_var()

    call init_case("10sp2", 5)
    allocate(Qs(n_hex))
    include "result_10sp2_5596079.inc"
    print *, tac(Qs)
    deallocate(Qs)
    call deallocate_var()

    print *, "sade ****************"
    np = 5*5*3*10
    rn(1) = 0.2d0
    call init_de("10sp1", 3, np)
    call init_population(np)
    print *, ymin
    print *, minloc(y_old), maxloc(y_old) 
    print *, minval(y_old), maxval(y_old) 
    call deallocate_de_var()
    call deallocate_var()

    print *, "sade for fixed topo ****************"
    idx_topo = (/8, 14, 20, 25, 42, 66, 76, 78/) ! 43392
    np = 5*5*4*10
    rn(1) = 0.2d0
    call init_de("10sp1", 4, np)
    call init_population_for_fixed_topo(n, idx_topo, np)
    print *, ymin
    print *, minloc(y_old), maxloc(y_old) 
    print *, minval(y_old), maxval(y_old) 
    print *, 8, 14, 20, 25, 42, 66, 76, 78
    do i=1, n_hex
        if(xmin(i,1)>0.1d-3) print *, i, xmin(i,1)
    enddo
    call deallocate_de_var()
    call deallocate_var()

    ! call run_sade("10sp1", 4, 5*5*4*10, 100, 100.d0, 100, 100, 0.3d0)

end program test 
