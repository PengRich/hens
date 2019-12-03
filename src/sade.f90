module sade 
    use utility
    use de_base
    implicit none
    integer(kind=4), parameter :: reset_std=1, debug=0
    integer(kind=4), parameter :: mutator_number=4, print_number=1000, &
        max_learning_period=1000, max_sampling_number=10000, &
        expected_number_err=5

    real(kind=8), public :: y_min_global, start, finish
    integer(kind=4), public :: n_hex_global, n_open, n_log_file

    real(kind=8), private, allocatable :: cf_normals(:, :), cr_normals(:, :), &
        cr_success(:, :), cf_success(:, :)
    real(kind=8), private :: crm(mutator_number), cfm(mutator_number), &
        crs(mutator_number), cfs(mutator_number), &
        prob_strategy(mutator_number), sum_prob_strategy
    integer(kind=4), private :: lp, n_lp, n_normal, sn, &
        success_strategy(mutator_number, max_learning_period), &
        failure_strategy(mutator_number, max_learning_period), &
        n_cfr(mutator_number), cfr_idx(mutator_number)
    contains
        subroutine get_strategy_id(id)
            implicit none
            real(kind=8) :: r, low, up 
            integer(kind=4) :: id
            r = rand(rn(1)) * sum_prob_strategy
            low = 0.d0
            up = 0.d0
            do id=1, mutator_number
                up = up + prob_strategy(id)
                if(r < up .and. r > low) exit
                low = low + prob_strategy(id)
            enddo
        end subroutine get_strategy_id

        subroutine update_prob_strategy()
            implicit none
            integer(kind=4) :: i
            real(kind=8) :: s(mutator_number), m, n, s0 

            do i=1, mutator_number
                m = sum(success_strategy(i, 1:lp))
                n = m+sum(failure_strategy(i, 1:lp))
                if(n>0) then
                    s(i) = real(m)/real(n)
                else
                    s(i) = 0.001d0
                endif
            enddo

            s0 = sum(s)
            do i=1, mutator_number
                prob_strategy(i) = max(0.1d0, s(i)/s0)
            enddo
            sum_prob_strategy = sum(prob_strategy)

        end subroutine update_prob_strategy

        ! subroutine update_crm(n)
        !     implicit none
        !     integer(kind=4) :: i, n
        !     real(kind=8) :: crs(n)

        !     do i=1, mutator_number
        !         crs = cr_success(i, :)
        !         call sort(n, crs)
        !         crm(i) = median(n, crs)
        !     enddo

        !     return
        ! end subroutine update_crm

        subroutine init_de_parameter_without_allocate()
            implicit none

            ! init parameter
            cf_normals = 0.3d0
            cr_normals = 0.5d0
            cr_success = 0.5d0
            cf_success = 0.5d0

            crm = 0.5d0
            crs = 0.3d0
            cfm = 0.5d0
            cfs = 0.3d0
            success_strategy = 0
            failure_strategy = 0
            prob_strategy = 1.d0 / real(mutator_number)
            sum_prob_strategy = sum(prob_strategy)
            n_lp = 0
            n_cfr = 1
            cfr_idx = 1
        end subroutine init_de_parameter_without_allocate

        subroutine init_de_parameter(np, sampling_number)
            implicit none
            integer(kind=4), intent(in) :: np, sampling_number

            n_normal = 2*np
            sn = min(sampling_number, max_sampling_number)

            ! allocate var
            allocate(cf_normals(mutator_number, n_normal))
            allocate(cr_normals(mutator_number, n_normal))
            allocate(cr_success(mutator_number, sn))
            allocate(cf_success(mutator_number, sn))

            ! init parameter
            cf_normals = 0.3d0
            cr_normals = 0.5d0
            cr_success = 0.5d0
            cf_success = 0.5d0

            crm = 0.5d0
            crs = 0.3d0
            cfm = 0.5d0
            cfs = 0.3d0
            success_strategy = 0
            failure_strategy = 0
            prob_strategy = 1.d0 / real(mutator_number)
            sum_prob_strategy = sum(prob_strategy)
            n_lp = 0
            n_cfr = 1
            cfr_idx = 1
        end subroutine init_de_parameter

        subroutine deallocate_de_parameter()
            implicit none
            deallocate(cf_normals)
            deallocate(cr_normals)
            deallocate(cf_success)
            deallocate(cr_success)
        end subroutine deallocate_de_parameter

        subroutine evolve(np, max_iter, qmin, learning_period)
            implicit none
            integer(kind=4), intent(in) :: np, max_iter, learning_period
            real(kind=8), intent(in) :: qmin

            integer(kind=4) :: i, j, k, m, strategy_id, &
                selected_strategy_record(np), expected_number
            real(kind=8) :: cf, cr, ep

            expected_number = sum(n_stms) + 1 + expected_number_err
            ! reset parameter
            lp = min(learning_period, max_learning_period)

            ! start iterate
            do i=1, max_iter
                u = x
                n_lp = n_lp + 1
                if(i > lp) then
                    call update_prob_strategy()
                    success_strategy(:, n_lp) = 0
                    failure_strategy(:, n_lp) = 0
                    do j=1, mutator_number
                        crm(j) = sum(cr_success(j, 1:lp)) / real(lp)
                        crs(j) = max(0.1d0, sqrt(sum((cr_success(j, 1:lp)- &
                            crm(j))**2.d0)/real(lp)))
                        cfm(j) = sum(cf_success(j, 1:lp)) / real(lp)
                        cfs(j) = max(0.1d0, sqrt(sum((cf_success(j, 1:lp)- &
                            cfm(j))**2.d0)/real(lp)))
                    enddo
                endif

                if(reset_std == 1 .and. mod(i, 10000) == 0) then
                    ! crm = 0.5d0
                    crs = 0.3d0
                    ! cfm = 0.5d0
                    cfs = 0.3d0
                endif
 
                do j=1, mutator_number
                    call generate_normal_rand(n_normal, rn(1), crm(j), crs(j), &
                        cr_normals(j, :))
                    call generate_normal_rand(n_normal, rn(1), cfm(j), cfs(j), &
                        cf_normals(j, :))
                enddo

                do j=1, np
                    call get_strategy_id(strategy_id)
                    selected_strategy_record(j) = strategy_id
                    do while(.true.)
                        if(cfr_idx(strategy_id) > n_normal) then
                            call generate_normal_rand(n_normal, rn(1), cfm(strategy_id), &
                                cfs(strategy_id), cf_normals(strategy_id, :)) 
                            call generate_normal_rand(n_normal, rn(1), crm(strategy_id), &
                                crs(strategy_id), cr_normals(strategy_id, :))
                            cfr_idx(strategy_id) = 1 
                        endif
                        cf = cf_normals(strategy_id, cfr_idx(strategy_id))
                        cr = cr_normals(strategy_id, cfr_idx(strategy_id))
                        cfr_idx(strategy_id) = cfr_idx(strategy_id) + 1
                        if(cr > 0.001d0 .and. cf > 0.001d0 .and. cr < 1.d0) exit 
                    enddo

                    ! cr = 0.8d0
                    ! cf = 0.5d0
                    ! strategy_id = 1
                    ! m = 1
                    do k=1, n_hex
                        select case(strategy_id)
                            case(2)
                                v(k, j) = random_to_best_two(np, j, k, cf)
                            case(3)
                                v(k, j) = random_one(np, j, k, cf) 
                            case(4)
                                v(k, j) = random_two(np, j, k, cf) 
                            case(5)
                                v(k, j) = current_to_random_one(np, j, k, cf) 
                            case(6)
                                v(k, j) = current_to_best_one(np, j, k, cf) 
                            !case(6)
                            !    v(k, j) = best_two(np, j, k, cf) 
                            !case(7)
                            !    v(k, j) = best_one(np, j, k, cf) 
                            case default
                                v(k, j) = best_to_random_one(np, j, k, cf) 
                        end select
                        ! if(v(k,j)<qmin .or. rand(rn(1)) < ep) v(k, j) = 0.d0
                        ! if(v(k,j)<qmin) v(k, j) = 0.d0
                        ! if(v(k,j)>0.1d-3) m = m+1
                    enddo

                    ! eliminate heat exchanger as number prob
                    ! ep = real(sum(n_stms)+4) / real(m)
                    ! if(m > sum(n_stms)+4) then
                    !     ep = real(min(elimination_number, m-sum(n_stms)-4)) / real(m)
                    ! else
                    !     ep = 1.d0 / real(m)
                    ! endif
                    ! ep = max(0.05d0, min(0.3d0, ep))
                    ! do k=1, n_hex
                    !     if(v(k,j)>0.1d-3 .and. rand(rn(1))<ep) v(k, j) = 0.d0
                    ! enddo

                    m = 0
                    do k=1, n_hex
                        if(rand(rn(1))<=cr .or. k==int(rand(rn(1))*real(n_hex))+1) u(k, j) = v(k, j)
                        if(u(k, j)>1.d-3) m = m + 1
                    enddo
                    if(m>expected_number) then
                        ep = real(expected_number) / real(m)
                    else
                        ! ep = real(m-1) / real(m)
                        ep = 1.d0
                    endif
                    do k=1, n_hex
                        if(u(k, j)>qmin) then
                            if(rand(rn(1))>ep) u(k,j) = 0.d0
                        else
                            u(k,j) = 0.d0
                        endif
                    enddo

                    y_new(j) = tac(u(:, j))
                ! enddo

                ! do j=1, np
                    strategy_id = selected_strategy_record(j)
                    if(y_new(j) < y_old(j)) then
                        y_old(j) = y_new(j)
                        x(:, j) = u(:, j)
                        success_strategy(strategy_id, n_lp) = success_strategy(strategy_id, n_lp) + 1
                        cr_success(strategy_id, n_cfr(strategy_id)) = cr 
                        cf_success(strategy_id, n_cfr(strategy_id)) = cf

                        n_cfr(strategy_id) = n_cfr(strategy_id) + 1
                        if(n_cfr(strategy_id) > sn) n_cfr(strategy_id) = 1
                    else
                        failure_strategy(strategy_id, n_lp) = failure_strategy(strategy_id, n_lp) + 1
                    endif
                enddo

                cfr_idx = 1 ! set normal random number index to 1

                ymin = minval(y_old)
                xmin = x(:, minloc(y_old))

                include 'print_sade.inc'

                if(n_lp >= lp) n_lp = 0

                if(sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np)) < 1.d0) then
                    exit
                endif

            enddo

            include 'print_sade_result.inc'
            return
        end subroutine evolve

        subroutine run_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number
            real(kind=8), intent(in) :: qmin !, elimination_prob 
            
            logical :: exist
            real(kind=8) :: random_state
            character(len=21) :: filename 
            n_open = 12
            n_log_file = 13
        
            filename = "output/" // trim(case_name) // "_sade.txt"
            call set_random_seed()

            do while(.true.)
                call cpu_time(start)
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(n_open, file=filename, status="old", position="append", action="write")
                else
                    open(n_open, file=filename, status="new", action="write")
                end if

                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1d0) exit
                enddo
                rn(1) = random_state
                call get_log_filename(case_name)
                open(unit=n_log_file, file=log_filename, action="write", status="replace")
                call init_de(case_name, stage, np)
                call init_population(np)
                call init_de_parameter(np, sampling_number)
                call evolve(np, max_iter, qmin, learning_period)
                call deallocate_de_parameter()
                call deallocate_de_var()
                call deallocate_var()
                call cpu_time(finish)
                close(n_log_file)
                write(*, *) log_filename, random_state, n_hex_global, y_min_global 
                write(n_open, *) log_filename, random_state, n_hex_global, &
                    y_min_global, finish-start
                close(n_open)
            enddo

            write(n_log_file, '("Time =", f10.1)') finish-start
            write(*, '("Time =", f10.1)') finish-start
            return
        end subroutine run_sade
end module sade 
