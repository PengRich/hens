module de_base
    use random_numbers
    use nosplit_simulator
    implicit none
    real(kind=8), public :: ymin
    real(kind=8), public, allocatable :: xmin(:, :)
    real(kind=8), public, allocatable :: x(:, :), v(:, :), u(:, :)
    real(kind=8), public, allocatable :: y_old(:), y_new(:) 
    real(kind=8), public :: rn(10)=0.1d0

    integer(kind=4), private :: max_match_number, i, j, k
    real(kind=8), private :: p_ex
    real(kind=8), private, allocatable :: ptnl(:, :), match_number(:, :)
    contains

        subroutine init_de(case_name, stage, np)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, np

            call init_case(case_name, stage)

            allocate(y_old(np))
            allocate(y_new(np))
            allocate(xmin(n_hex, 1))
            allocate(x(n_hex, np))
            allocate(v(n_hex, np))
            allocate(u(n_hex, np))
            y_old = 1.d20
            y_new = 1.d20
            xmin = 0.d0
            ymin = 1.d20
            x = 0.d0
            u = 0.d0
            v = 0.d0
        end subroutine init_de

        subroutine init_population(np)
            implicit none
            integer(kind=4), intent(in) :: np
            integer(kind=4) :: hs, cs, st, max_try
            ! real(kind=8), private :: rn(3)=0.1d0

            max_match_number = int(maxval(n_stms)+1)
            p_ex = real(sum(n_stms)+2) / real(n_hex)
            allocate(ptnl(2, maxval(n_stms)))
            allocate(match_number(2, maxval(n_stms)))
            do i=1, np
                ! x(:, i) = 0.d0
                match_number = 0
                ptnl = stms%ptnl

                k = 0
                do j=1, n_hex 
                    if(rand(rn(1))<p_ex) k = k + 1
                enddo

                max_try = 0
                do while(k > 0)
                    hs = int(rand(rn(1))*n_stms(1)+1)
                    cs = int(rand(rn(1))*n_stms(2)+1)
                    st = int(rand(rn(1))*n_st+1)
                    j = unit_id(hs, cs, st)
                    if(abs(x(j, i)) > 1.d-3) cycle
                    if(max_try < 1000) then
                        if(match_number(1, idx(j)%hstm_id) > max_match_number) cycle
                        if(match_number(2, idx(j)%cstm_id) > max_match_number) cycle
                    ! else:
                    !    max_try = 0
                    endif

                    match_number(1, idx(j)%hstm_id) = match_number(1, idx(j)%hstm_id) + 1
                    match_number(2, idx(j)%cstm_id) = match_number(2, idx(j)%cstm_id) + 1

                    x(j, i) = max(1.d0, rand(rn(1))*(min(ptnl(1, idx(j)%hstm_id), ptnl(2, idx(j)%cstm_id))))
                    ptnl(1, idx(j)%hstm_id) = ptnl(1, idx(j)%hstm_id) - x(j, i)
                    ptnl(2, idx(j)%cstm_id) = ptnl(2, idx(j)%cstm_id) - x(j, i)

                    k = k - 1
                    max_try = max_try + 1
                enddo
                y_old(i) = tac(x(:, i))
            enddo

            ymin = minval(y_old)
            xmin = x(:, minloc(y_old))
 
            deallocate(ptnl)
            deallocate(match_number)
        end subroutine init_population

        subroutine deallocate_de_var
            implicit none
            deallocate(y_old)
            deallocate(y_new)
            deallocate(xmin)
            deallocate(x)
            deallocate(v)
            deallocate(u)
        end subroutine deallocate_de_var

        subroutine init_population_for_fixed_topo(n_topo, idx_topo, np)
            implicit none
            integer(kind=4) :: n_topo, idx_topo(n_topo), np

            allocate(ptnl(2, int(maxval(n_stms))))

            do i=1, np
                x(:, i) = 0.d0
                ptnl = stms%ptnl

                do k=1, n_topo 
                    j = idx_topo(k) 
                    x(j, i) = max(1.d0, rand(rn(1))*(min(ptnl(1, idx(j)%hstm_id), ptnl(2, idx(j)%cstm_id))))
                    ptnl(1, idx(j)%hstm_id) = ptnl(1, idx(j)%hstm_id) - x(j, i)
                    ptnl(2, idx(j)%cstm_id) = ptnl(2, idx(j)%cstm_id) - x(j, i)

                enddo
                y_old(i) = tac(x(:, i))
            enddo

            ymin = minval(y_old)
            xmin = x(:, minloc(y_old))
 
            deallocate(ptnl)
        end subroutine init_population_for_fixed_topo

        function random_one(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(3)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            if(ids(1) == j) then
                ids(1) = int(rand(rn(1))*np) + 1
            endif
            ids(2) = int(rand(rn(1))*np) + 1
            if(ids(2)==ids(1) .or. ids(2)==j) then
                ids(2) = int(rand(rn(1))*np) + 1
            endif
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
            vec = x(k, ids(1)) + cf*(x(k, ids(2))-x(k, ids(3)))

        end function random_one 

        function best_one(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(2)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1)*np)) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
 
            vec = xmin(k, 1) + cf*(x(k, ids(1))-x(k, ids(2)))

        end function best_one

        function current_to_best_one(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(2)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1)*np)) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
 
            vec = x(k, j) + cf*(xmin(k, 1)-x(k, j)) + cf*(x(k, ids(1))-x(k, ids(2)))

        end function current_to_best_one 

        function current_to_random_one(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(3)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1)*np)) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
 
            vec = x(k, j) + rand(rn(1))*(x(k, ids(1))-x(k, j)) + cf*(x(k, ids(2))-x(k, ids(3)))

        end function current_to_random_one 

        function best_two(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(4)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1)*np)) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
            ids(4) = int(rand(rn(1))*np) + 1
            do while(ids(4)==ids(1) .or. ids(4)==ids(2) .or. ids(4)==ids(3) .or.ids(4)==j)
                ids(4) = int(rand(rn(1))*np) + 1
            enddo

            vec = xmin(k, 1) + cf*(x(k, ids(1))-x(k, ids(2))) + cf*(x(k, ids(3))-x(k, ids(4)))

        end function best_two

        function random_two(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(5)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            if(ids(1) == j) then
                ids(1) = int(rand(rn(1))*np) + 1
            endif
            ids(2) = int(rand(rn(1))*np) + 1
            if(ids(2)==ids(1) .or. ids(2)==j) then
                ids(2) = int(rand(rn(1))*np) + 1
            endif
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
            ids(4) = int(rand(rn(1))*np) + 1
            do while(ids(4)==ids(1) .or. ids(4)==ids(2) .or. ids(4)==ids(3) .or.ids(4)==j)
                ids(4) = int(rand(rn(1))*np) + 1
            enddo
            ids(5) = int(rand(rn(1))*np) + 1
            do while(ids(5)==ids(1) .or. ids(5)==ids(2) .or. ids(5)==ids(3) .or. ids(5)==ids(4) .or.ids(5)==j)
                ids(5) = int(rand(rn(1))*np) + 1
            enddo

            vec = x(k, ids(1)) + cf*(x(k, ids(2))-x(k, ids(3))) + cf*(x(k, ids(4))-x(k, ids(5)))

        end function random_two 

        function random_to_best_two(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(4)
            real(kind=8) :: vec

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1)*np)) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
            ids(4) = int(rand(rn(1))*np) + 1
            do while(ids(4)==ids(1) .or. ids(4)==ids(2) .or. ids(4)==ids(3) .or.ids(4)==j)
                ids(4) = int(rand(rn(1))*np) + 1
            enddo
 
            vec = x(k, j) + cf*(xmin(k, 1)-x(k, j)) + cf*(x(k, ids(1))-x(k, ids(2))) &
                + cf*(x(k, ids(3))-x(k, ids(4)))

        end function random_to_best_two

        function best_to_random_one(np, j, k, cf) result(vec)
            implicit none
            integer(kind=4), intent(in) :: np, j, k
            real(kind=8), intent(in) :: cf
            integer(kind=4) :: ids(4)
            real(kind=8) :: vec, cf0

            ids(1) = int(rand(rn(1))*np) + 1
            do while(ids(1) == j)
                ids(1) = int(rand(rn(1))*np) + 1
            enddo
            ids(2) = int(rand(rn(1))*np) + 1
            do while(ids(2)==ids(1) .or. ids(2)==j)
                ids(2) = int(rand(rn(1))*np) + 1
            enddo
            ids(3) = int(rand(rn(1))*np) + 1
            do while(ids(3)==ids(1) .or. ids(3)==ids(2) .or. ids(3)==j)
                ids(3) = int(rand(rn(1))*np) + 1
            enddo
            ids(4) = int(rand(rn(1))*np) + 1
            do while(ids(4)==ids(1) .or. ids(4)==ids(2) .or. ids(4)==ids(3) .or.ids(4)==j)
                ids(4) = int(rand(rn(1))*np) + 1
            enddo
 
            cf0 = rand(rn(1))
            vec = xmin(k, 1) + cf*(x(k, ids(1))-x(k, ids(2))) + &
                cf0*(x(k, ids(3))-x(k, ids(4)))

        end function best_to_random_one 

end module de_base
