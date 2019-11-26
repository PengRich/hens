! 1-hot, 2-cold
module types 
    implicit none 
    type :: stream
        real(kind=8) :: T_src, T_trg, F, h, ptnl
    end type stream

    type :: cost_factor 
        real(kind=8) :: hu, cu, f, m, n
    end type cost_factor 

    type :: heat_exchanger
        integer(kind=4) :: sta
        real(kind=8) :: A, htc, cost
    end type heat_exchanger

    type :: semi_heat_exchanger
        integer(kind=4) :: id
        real(kind=8) :: T(2) ! 1-in, 2-out
        real(kind=8) :: A, htcs(2), ucfs(2)
        type(stream), pointer :: stm, ustms(:)
        type(heat_exchanger), pointer :: parent
        type(semi_heat_exchanger), pointer :: match, next
    end type semi_heat_exchanger

    type :: utility_heat_exchanger 
        real(kind=8) :: T, Q, A, htcs(2), ucfs(2), cost
    end type utility_heat_exchanger 

    type :: heat_exchanger_index
        integer(kind=4) :: id, hstm_id, cstm_id, stage_id
        type(heat_exchanger_index), pointer :: next
    end type heat_exchanger_index 

end module types

module nosplit_simulator
    use types
    implicit none
    !***************************
    ! Public
    !***************************
    integer(kind=4), public :: n_stms(2), n_st, n_hex
    type(stream), public, target :: ustms(2)
    type(heat_exchanger_index), public, allocatable, target :: idx(:)
    type(heat_exchanger_index), public, allocatable, target :: idx0(:, :, :)
    type(stream), public, allocatable, target :: stms(:, :)
    !***************************
    ! Private 
    !***************************
    type(cost_factor), private :: factor
    type(heat_exchanger), private, allocatable, target :: hexs(:)
    type(semi_heat_exchanger), private, allocatable, target :: heads(:, :) 
    type(semi_heat_exchanger), private, allocatable, target :: shexs(:, :, :, :)
    type(utility_heat_exchanger), private, allocatable, target :: uhexs(:, :)
    !***************************
    ! Private, function 
    !***************************
    private :: heat_transfer_area
    private :: penalty 
    private :: heat_transfer_coefficent 
    private :: init_topo
    private :: cal_grid_temp 
    private :: cal_internal_cost 
    private :: cal_terminal_cost 
    contains
        !***************************
        ! Calculate index
        !***************************
        function unit_id(i, j, k) result(id)
            implicit none
            integer(kind=4), intent(in) :: i, j, k 
            integer(kind=4) :: id

            id = (k-1)*product(n_stms) + (i-1)*n_stms(2) + j 

            return
        end function unit_id 
        !***************************
        ! Heat tranfer coefficient 
        !***************************
        function heat_transfer_coefficent(h1, h2) result(htc)
            implicit none
            real(kind=8), intent(in) :: h1, h2
            real(kind=8) :: htc

            htc = (h1*h2) / (h1+h2)

            return
        end function heat_transfer_coefficent
        !***************************
        ! Heat tranfer area 
        !***************************
        function heat_transfer_area(Q, htc, dtl, dtr) result(A)
            implicit none
            real(kind=8), intent(in) :: Q, htc, dtl, dtr
            real(kind=8) :: dtm, A
            
            if(abs(dtl-dtr) .gt. 1.d-5) then
                dtm = abs(dtl-dtr)/abs(log(dtl/dtr))
            else
                dtm = abs(dtl+dtr)/2.d0
            endif

            A = abs(Q)/dtm/htc

            return
        end function heat_transfer_area 
        !***************************
        ! Penalty function 
        !***************************
        real(kind=8) function penalty(dev)
            implicit none
            real(kind=8), intent(in) :: dev
            real(kind=8) :: penalty_factor

            penalty_factor = 1.d7
            penalty = 0.5d0*penalty_factor*(max(0.d0, -dev+1.d0)**2.d0)

            return 
        end function penalty
        !***************************
        ! Allocate
        !***************************
        subroutine allocate_var()
            implicit none
            n_hex = product(n_stms) * n_st
            allocate(idx(n_hex))
            allocate(idx0(n_stms(1), n_stms(2), n_st))
            allocate(stms(2, maxval(n_stms)))
            allocate(heads(2, maxval(n_stms)))
            allocate(hexs(n_hex))
            allocate(shexs(2, n_stms(1), n_stms(2), n_st))
            allocate(uhexs(2, maxval(n_stms)))
        end subroutine allocate_var

        subroutine deallocate_var()
            implicit none
            deallocate(idx)
            deallocate(idx0)
            deallocate(stms)
            deallocate(heads)
            deallocate(hexs)
            deallocate(shexs)
            deallocate(uhexs)
        end subroutine deallocate_var
        !***************************
        ! Initialize case parameter 
        !***************************
        subroutine init_case(case_name, stage)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage 
            integer(kind=4) :: i, j, k, id
            real(kind=8) :: f

            n_st = stage
            select case(case_name)
                case("4sp1")
                    n_stms = (/2, 2/)
                    call allocate_var()
                    include "4sp1.inc"
                case("9sp1")
                    n_stms = (/4, 5/)
                    call allocate_var()
                    include "9sp1.inc"
                case("10sp1")
                    n_stms = (/5, 5/)
                    call allocate_var()
                    include "10sp1.inc"
                case("10sp2")
                    n_stms = (/6, 4/) 
                    call allocate_var()
                    include "10sp2.inc"
                case("15sp1")
                    n_stms = (/8, 7/)
                    call allocate_var()
                    include "15sp1.inc"
                case("20sp1")
                    n_stms = (/10, 10/)
                    call allocate_var()
                    include "20sp1.inc"
                case default
                    stop "Input wrong case name"
            end select

            do k=1, n_st
                do i=1, n_stms(1)
                    do j=1, n_stms(2)
                        id = unit_id(i, j, k)
                        shexs(1, i, j, k)%id = id
                        shexs(1, i, j, k)%stm => stms(1, i)
                        shexs(1, i, j, k)%htcs(1) = heat_transfer_coefficent(stms(1, i)%h, ustms(2)%h)
                        shexs(1, i, j, k)%htcs(2) = heat_transfer_coefficent(stms(1, i)%h, ustms(1)%h)
                        shexs(1, i, j, k)%ucfs(1) = factor%cu
                        shexs(1, i, j, k)%ucfs(2) = factor%hu
                        shexs(1, i, j, k)%ustms => ustms(2:1:-1) 
                        shexs(1, i, j, k)%parent => hexs(id)
                        shexs(1, i, j, k)%match => shexs(2, i, j, k)

                        shexs(2, i, j, k)%id = id
                        shexs(2, i, j, k)%stm => stms(2, j)
                        shexs(2, i, j, k)%htcs(1) = heat_transfer_coefficent(stms(2, j)%h, ustms(1)%h)
                        shexs(2, i, j, k)%htcs(2) = heat_transfer_coefficent(stms(2, j)%h, ustms(2)%h)
                        shexs(2, i, j, k)%ucfs(1) = factor%hu
                        shexs(2, i, j, k)%ucfs(2) = factor%cu
                        shexs(2, i, j, k)%ustms => ustms
                        shexs(2, i, j, k)%parent => hexs(id) 
                        shexs(2, i, j, k)%match => shexs(1, i, j, k) 

                        hexs(id)%htc = heat_transfer_coefficent(stms(1, i)%h, stms(2, j)%h)

                        idx(id)%id = id
                        idx(id)%hstm_id = i 
                        idx(id)%cstm_id = j 
                        idx0(i, j, k)%id = id
                        idx0(i, j, k)%hstm_id = i 
                        idx0(i, j, k)%cstm_id = j 

                    enddo
                enddo
            enddo

            do k=1, 2
                do i=1, n_stms(k)
                    uhexs(k, i)%htcs(1) = heat_transfer_coefficent(stms(k, i)%h, ustms(1)%h)
                    uhexs(k, i)%htcs(2) = heat_transfer_coefficent(stms(k, i)%h, ustms(2)%h)
                    uhexs(k, i)%ucfs(1) = factor%hu
                    uhexs(k, i)%ucfs(2) = factor%cu
                    stms(k, i)%ptnl = abs(stms(k, i)%F*(stms(k, i)%T_src-stms(k, i)%T_trg))
                enddo
            enddo

            do i=1, n_stms(1) 
                if(stms(1, i)%T_trg < ustms(2)%T_trg) then
                    stms(1, i)%ptnl = stms(1, i)%F*(stms(1, i)%T_src-(ustms(2)%T_trg+1.d0))
                endif
            enddo
            
            do j=1, n_stms(2) 
                if(stms(2, j)%T_trg > ustms(1)%T_trg) then
                    stms(2, j)%ptnl = stms(2, j)%F*((ustms(1)%T_trg-1.d0)-stms(2, j)%T_src)
                endif
            enddo

            return
        end subroutine init_case
        !***************************
        ! Figure out topology 
        !***************************
        subroutine init_topo(Qs)
            implicit none
            real(kind=8), intent(in) :: Qs(n_hex)
            integer(kind=4) :: i, j, k
            type(semi_heat_exchanger), pointer:: next 

            do i=1, n_stms(1) 
                next => heads(1, i)
                do k=1, n_st
                    do j=1, n_stms(2)
                        if(abs(Qs(idx0(i, j, k)%id)) < 1.d-3) cycle
                        next%next => shexs(1, i, j, k)
                        next => shexs(1, i, j, k)
                    enddo
                enddo
                next%next => null()
            enddo

            do j=1, n_stms(2) 
                next => heads(2, j)
                do k=n_st, 1, -1
                    do i=n_stms(1), 1, -1
                        if(abs(Qs(idx0(i, j, k)%id)) < 1.d-3) cycle
                        next%next => shexs(2, i, j, k)
                        next => shexs(2, i, j, k)
                    enddo
                enddo
                next%next => null()
            enddo

            return
        end subroutine init_topo
        !***************************
        ! Calculate grid temperature
        !***************************
        subroutine cal_grid_temp(Qs)
            implicit none
            real(kind=8), intent(in) :: Qs(n_hex)
            integer(kind=4) :: i, j, k
            real(kind=8) :: f
            logical :: cond(2)
            type(semi_heat_exchanger), pointer :: next 

            uhexs%Q = 0.d0
            f = 1.d0
            do k=1, 2
                uhexs(k, :)%T = stms(k, :)%T_src
                do i=1, n_stms(k) 
                    if(associated(heads(k, i)%next)) then
                        next => heads(k, i)%next
                        do while(.true.)
                            next%T(1) = uhexs(k, i)%T
                            next%T(2) = next%T(1) - f*(Qs(next%id)/stms(k, i)%F)
                            uhexs(k, i)%T = next%T(2)
                            if(.not. associated(next%next)) exit 
                            next => next%next
                        enddo
                    endif
                    uhexs(k, i)%Q = (uhexs(k, i)%T-stms(k, i)%T_trg) * stms(k, i)%F
                enddo
                f = -1.d0
            enddo

            do i=1, n_stms(1)
                if(.not. associated(heads(1, i)%next)) cycle
                next => heads(1, i)%next
                do while(.true.)
                    cond(1) = next%T(1) < next%match%T(2)
                    cond(2) = next%T(2) < next%match%T(1)
                    next%parent%sta = 3
                    if(Qs(next%id)>0.d0) then
                        if(cond(1) .or. cond(2)) next%parent%sta = 1
                    else
                        if(.not.cond(1) .or. .not.cond(2)) next%parent%sta = 2
                    endif
                    if(.not. associated(next%next)) exit
                    next => next%next
                enddo
            enddo

            return
        end subroutine cal_grid_temp
        !***************************
        ! Calculate internal cost
        !***************************
        subroutine cal_internal_cost(Qs)
            implicit none
            real(kind=8), intent(in) :: Qs(n_hex)
            integer(kind=4) :: i, sta
            real(kind=8) :: dtl, dtr, f(2)
            type(semi_heat_exchanger), pointer:: next 

            f = (/1.d0, -1.d0/)

            hexs%cost = 0.d0
            do i=1, n_stms(1) 
                if(.not. associated(heads(1, i)%next)) cycle
                next => heads(1, i)%next
                do while(.true.)
                    sta = next%parent%sta
                    if(sta==3) then
                        dtl = next%T(1) - next%match%T(2)
                        dtr = next%T(2) - next%match%T(1)
                        next%parent%A = heat_transfer_area(Qs(next%id), next%parent%htc, dtl, dtr)
                        next%parent%cost = factor%f + factor%m*(next%parent%A**factor%n)
                    else
                        dtl = f(sta) * (next%T(1)-next%ustms(sta)%T_trg)
                        dtr = f(sta) * (next%T(2)-next%ustms(sta)%T_src)
                        if(dtl>0.d0 .and. dtr>0.d0) then
                            next%A = heat_transfer_area(Qs(next%id), next%htcs(sta), dtl, dtr)
                            next%parent%cost = factor%f + factor%m*(next%A**factor%n)
                            next%parent%cost = next%parent%cost + next%ucfs(sta)*abs(Qs(next%id))
                        else
                            next%parent%cost = penalty(dtl) + penalty(dtr)
                        endif
                        dtl = f(sta) * (next%match%ustms(sta)%T_trg-next%match%T(1))
                        dtr = f(sta) * (next%match%ustms(sta)%T_src-next%match%T(2))
                        if(dtl>0.d0 .and. dtr>0.d0) then
                            next%match%A = heat_transfer_area(Qs(next%id), next%match%htcs(sta), dtl, dtr)
                            next%parent%cost = next%parent%cost + factor%f + factor%m*(next%match%A**factor%n)
                            next%parent%cost = next%parent%cost + next%match%ucfs(sta)*abs(Qs(next%id))
                        else
                            next%parent%cost = next%parent%cost + penalty(dtl) + penalty(dtr)
                        endif
                    endif
                    if(.not. associated(next%next)) exit
                    next => next%next
                enddo
            enddo

            return
        end subroutine cal_internal_cost 
        !***************************
        ! Calculate terminal cost
        !***************************
        subroutine cal_terminal_cost()
            implicit none
            integer(kind=4) :: sta, i, k
            real(kind=8) :: dtl, dtr, f

            uhexs%cost = 0.d0
            do k=1, 2
                do i=1, n_stms(k)
                    if(abs(uhexs(k, i)%Q)<1.d-3) cycle
                    if(uhexs(k, i)%Q>0.d0) then
                        sta = 2
                        f = 1.d0
                    else
                        sta = 1 
                        f = -1.d0
                    endif
                    dtl = f*(uhexs(k, i)%T - ustms(sta)%T_trg)
                    dtr = f*(stms(k, i)%T_trg - ustms(sta)%T_src)
                    if(dtl>0.d0 .and. dtr>0.d0) then
                        uhexs(k, i)%A = heat_transfer_area(uhexs(k, i)%Q, uhexs(k, i)%htcs(sta), dtl, dtr) 
                        uhexs(k, i)%cost = factor%f + factor%m*(uhexs(k, i)%A**factor%n)
                        uhexs(k, i)%cost = uhexs(k, i)%cost + uhexs(k, i)%ucfs(sta)*abs(uhexs(k, i)%Q)
                    else
                        uhexs(k, i)%cost = penalty(dtl) + penalty(dtr)
                    endif
                enddo
            enddo

            return
        end subroutine cal_terminal_cost

        real(kind=8) function tac0(Qs)
            implicit none
            real(kind=8) :: Qs(n_hex)
  
            call cal_grid_temp(Qs)
            call cal_internal_cost(Qs) 
            call cal_terminal_cost()
            tac0 = sum(hexs%cost) + sum(uhexs%cost) 

            return
        end function tac0 

        real(kind=8) function tac(Qs)
            implicit none
            real(kind=8), intent(in) :: Qs(n_hex)

            call init_topo(Qs)
            tac = tac0(Qs)
            return
        end function tac 

        subroutine print_result(Qs)
            implicit none
            real(kind=8), intent(in) :: Qs(n_hex)
            integer(kind=4) :: i
            write(*, *) "result"
            do i=1, n_hex
                if(abs(Qs(i))>1.d-3) write(*, *) "Q[", i, "]=", Qs(i)
            enddo
            write(*, *) "tac=", tac(Qs)
        end subroutine print_result

end module nosplit_simulator
