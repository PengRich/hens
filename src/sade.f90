module sade 
    use utility
    use de_base
    implicit none
    integer(kind=4), parameter :: reset_std=1, debug=0
    integer(kind=4), parameter :: mutator_number=4, print_number=1000, &
        max_learning_period=1000, max_sampling_number=10000

    real(kind=8), private, allocatable :: cf_normals(:, :), cr_normals(:, :), &
        cr_success(:, :), cf_success(:, :)
    real(kind=8), private :: start, finish, &
        crm(mutator_number), cfm(mutator_number), &
        crs(mutator_number), cfs(mutator_number), &
        prob_strategy(mutator_number), sum_prob_strategy, &
        y_min_global
    integer(kind=4), private :: lp, n_lp, n_hex_global, &
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
                if(id==5) then
                    print *, sum_prob_strategy, low-prob_strategy(4), r, up
                    stop
                endif


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

        subroutine evolve(np, max_iter, qmin, learning_period, sampling_number, &
                elimination_prob)
            implicit none
            integer(kind=4), intent(in) :: np, max_iter, learning_period, &
                sampling_number
            real(kind=8), intent(in) :: qmin, elimination_prob

            integer(kind=4) :: i, j, k, sn, strategy_id, &
                selected_strategy_record(np), n_normal
            real(kind=8) :: cf, cr, ep

            n_normal = 2*np

            call cpu_time(start)
            ! reset parameter
            lp = min(learning_period, max_learning_period)
            sn = min(sampling_number, max_sampling_number)
            ep = min(1.d0, max(0.d0, elimination_prob))

            ! allocate var
            allocate(cf_normals(mutator_number, n_normal))
            allocate(cr_normals(mutator_number, n_normal))
            allocate(cr_success(mutator_number, sn))
            allocate(cf_success(mutator_number, sn))

            ! init parameter
            cf_normals = 0.3d0
            cr_normals = 0.5d0

            crm = 0.5d0
            crs = 0.3d0
            cfm = 0.5d0
            cfs = 0.3d0
            success_strategy = 0
            failure_strategy = 0
            cr_success = 0.5d0
            cf_success = 0.5d0
            prob_strategy = 1.d0 / real(mutator_number)
            sum_prob_strategy = sum(prob_strategy)
            n_lp = 0
            n_cfr = 1
            cfr_idx = 1

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
                            print *, mutator_number, strategy_id, cfr_idx(strategy_id), "here"
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
                        if(v(k,j)<qmin .or. rand(rn(1)) < ep) v(k, j) = 0.d0
                    enddo

                    do k=1, n_hex
                        if(rand(rn(1))<=cr .or. k==int(rand(rn(1))*real(n_hex))+1) u(k, j) = v(k, j)
                    enddo
                    y_new(j) = tac(u(:, j))
                enddo

                do j=1, np
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

            call cpu_time(finish)
            include 'print_sade_result.inc'
            deallocate(cf_normals)
            deallocate(cr_normals)
            deallocate(cf_success)
            deallocate(cr_success)
            return
        end subroutine evolve

        subroutine run_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, elimination_prob)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number
            real(kind=8), intent(in) :: qmin, elimination_prob 
            
            logical :: exist
            real(kind=8) :: random_state
            character(len=21) :: filename 
            integer :: values(1:8), k, i, j
            integer, allocatable :: seed(:)
        
            filename = "output/" // trim(case_name) // "_sade.txt"
            call date_and_time(values=values)
            call random_seed(size=k)
            allocate(seed(1:k))
            seed(:) = values(8)
            call random_seed(put=seed)

            do while(.true.)
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(12, file=filename, status="old", position="append", action="write")
                else
                    open(12, file=filename, status="new", action="write")
                end if

                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1d0) exit
                enddo
                rn(1) = random_state
                call get_log_filename(case_name)
                open(unit=13, file=log_filename, action="write", status="replace")
                call init_de(case_name, stage, np)
                call init_population(np)
                call evolve(np, max_iter, qmin, learning_period, &
                    sampling_number, elimination_prob)
                call deallocate_de_var()
                call deallocate_var()
                write(*, *) log_filename, random_state, n_hex_global, y_min_global 
                write(12, *) log_filename, random_state, n_hex_global, y_min_global 
                close(12)
                close(13)
                deallocate(log_filename)
            enddo

            deallocate(seed)
            return
        end subroutine run_sade

        ! subroutine evolve_fixed(n, idx_q, np, max_iter, qmin, lp0, prob)
        !     implicit none
        !     integer(kind=4) :: n, idx_q(n)
        !     integer(kind=4) :: np, max_iter, lp0
        !     real(kind=8) :: qmin, prob

        !     integer(kind=4) :: n_lp, n_cfr(mutator_number), cfr_idx(mutator_number)
        !     integer(kind=4) :: i, j, k, strategy_id, lp, k0
        !     real(kind=8) :: cr_success(mutator_number, lp0), cf_success(mutator_number, lp0)
        !     real(kind=8) :: cf, cr
        !     real(kind=8) :: crm(mutator_number), cfm(mutator_number)
        !     real(kind=8) :: crs(mutator_number), cfs(mutator_number)
        !     real(kind=8), allocatable :: cf_normals(:, :), cr_normals(:, :)

        !     call cpu_time(start)
        !     call get_log_filename()
        !     open(unit=20, file=log_filename, action="write", status="replace")
        !     ! allocate var
        !     allocate(cf_normals(mutator_number, np))
        !     allocate(cr_normals(mutator_number, np))

        !     cf_normals = 0.3d0
        !     cr_normals = 0.5d0

        !     lp = min(lp0, max_learning_period)
        !     ! init parameter
        !     crm = 0.5d0
        !     crs = 0.3d0
        !     cfm = 0.5d0
        !     cfs = 0.3d0
        !     success_strategy = 0
        !     failure_strategy = 0
        !     cr_success = 0.5d0
        !     cf_success = 0.5d0
        !     prob_strategy = 1.d0 / real(mutator_number)
        !     n_lp = 0
        !     n_cfr = 1
        !     cfr_idx = 1

        !     do i=1, max_iter
        !         n_lp = n_lp + 1
        !         if(i > lp) then
        !             call update_prob_strategy(i, lp)
        !             success_strategy(:, n_lp) = 0
        !             failure_strategy(:, n_lp) = 0
        !             crm(strategy_id) = sum(cr_success(strategy_id, 1:lp)) / real(lp)
        !             crs(strategy_id) = max(0.1d0, sqrt(sum((cr_success(strategy_id, 1:lp)- &
        !                 crm(strategy_id))**2.d0)/real(lp)))
        !             cfm(strategy_id) = sum(cf_success(strategy_id, 1:lp)) / real(lp)
        !             cfs(strategy_id) = max(0.1d0, sqrt(sum((cf_success(strategy_id, 1:lp)- &
        !                 cfm(strategy_id))**2.d0)/real(lp)))
        !         endif

        !         if(reset_std == 1 .and. mod(i, 1000) == 0) then
        !             ! crm = 0.5d0
        !             crs = 0.3d0
        !             ! cfm = 0.5d0
        !             cfs = 0.3d0
        !         endif
 
        !         do j=1, mutator_number
        !             call generate_normal_rand(np, rn(1), cfm(j), cfs(j), cf_normals(j, :)) 
        !             call generate_normal_rand(np, rn(1), crm(j), crs(j), cr_normals(j, :))
        !         enddo

        !         do j=1, np
        !             call get_strategy_id(strategy_id)
        !             do while(.true.)
        !                 if(cfr_idx(strategy_id) == np+1) then
        !                     call generate_normal_rand(np, rn(1), cfm(strategy_id), cfs(strategy_id), cf_normals(strategy_id, :)) 
        !                     call generate_normal_rand(np, rn(1), crm(strategy_id), crs(strategy_id), cr_normals(strategy_id, :))
        !                     cfr_idx(strategy_id) = 1 
        !                 endif
        !                 cf = cf_normals(strategy_id, cfr_idx(strategy_id))
        !                 cr = cr_normals(strategy_id, cfr_idx(strategy_id))
        !                 cfr_idx(strategy_id) = cfr_idx(strategy_id) + 1
        !                 if(cr > 0.001d0 .and. cf > 0.001d0 .and. cr < 1.d0) exit 
        !             enddo

        !             ! cr = 0.8d0
        !             ! cf = 0.5d0
        !             ! strategy_id = 5
        !             do k0=1, n 
        !                 k = idx_q(k0)
        !                 select case(strategy_id)
        !                     case(2)
        !                         v(k, j) = random_to_best_two(np, j, k, cf)
        !                     case(3)
        !                         v(k, j) = random_one(np, j, k, cf) 
        !                     case(4)
        !                         v(k, j) = random_two(np, j, k, cf) 
        !                     case(5)
        !                         v(k, j) = current_to_random_one(np, j, k, cf) 
        !                     case(6)
        !                         v(k, j) = current_to_best_one(np, j, k, cf) 
        !                     case(7)
        !                         v(k, j) = best_two(np, j, k, cf) 
        !                     case(8)
        !                         v(k, j) = best_one(np, j, k, cf) 
        !                     case default
        !                         v(k, j) = best_to_random_one(np, j, k, cf) 
        !                 end select
        !                 if(v(k,j)<qmin .or. rand(rn(1)) < prob) v(k, j) = 0.d0
        !             enddo

        !             do k=1, n_hex
        !                 if(rand(rn(1))<=cr .or. k==int(rand(rn(1))*n_hex)+1) u(k, j) = v(k, j)
        !             enddo

        !             y_new(j) = tac(u(:, j))

        !             if(y_new(j) < y_old(j)) then
        !                 y_old(j) = y_new(j)
        !                 x(:, j) = u(:, j)

        !                 success_strategy(strategy_id, n_lp) = success_strategy(strategy_id, n_lp) + 1
        !                 cr_success(strategy_id, n_cfr(strategy_id)) = cr 
        !                 cf_success(strategy_id, n_cfr(strategy_id)) = cf

        !                 n_cfr(strategy_id) = n_cfr(strategy_id) + 1
        !                 if(n_cfr(strategy_id) == lp+1) n_cfr(strategy_id) = 1
        !             else
        !                 failure_strategy(strategy_id, n_lp) = failure_strategy(strategy_id, n_lp) + 1
        !             endif
        !         enddo

        !         cfr_idx = 1
        !         ymin = minval(y_old)
        !         xmin = x(:, minloc(y_old))
        !         
        !         if(mod(i, print_number) == 0) then
        !             k = 0
        !             do j=1, n_hex
        !                 if(xmin(j, 1)>1.d-3) k = k+1
        !             enddo
        !             if(debug==1) then
        !                 write(20, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        !                 write(20, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        !                 write(20, *) i, "p", prob_strategy, sum(prob_strategy)
        !                 write(20, *) i, "crm", crm
        !                 write(20, *) i, "crs", crs
        !                 write(20, *) i, "cfm", cfm
        !                 write(20, *) i, "cfs", cfs
        !                 write(*, *) i, ymin, k, sum(y_old)/real(np), maxval(y_old)
        !                 write(*, *) i, minloc(y_old), maxloc(y_old), sum(success_strategy(:, n_lp))
        !                 write(*, *) i, "p", prob_strategy, sum(prob_strategy)
        !                 write(*, *) i, "crm", crm
        !                 write(*, *) i, "crs", crs
        !                 write(*, *) i, "cfm", cfm
        !                 write(*, *) i, "cfs", cfs
        !                 write(*, *) i, "n_cfr", n_cfr 
        !             else
        !                 write(20, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
        !                     sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
        !                 write(*, "(1x, i10, i5, 2x, f15.1, 2x, f15.1, 2x, f15.1, 2x, f15.1)") i, k, ymin, &
        !                     sum(y_old)/real(np), maxval(y_old), sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np))
        !             endif
        !         endif
        !         if(n_lp == lp) n_lp = 0
        !         if(sqrt(sum((y_old - sum(y_old)/real(np))**2.d0)/real(np)) < 1.d0) then
        !             exit
        !         endif
        !     enddo

        !     do i=1, n_hex
        !         if(abs(xmin(i, 1))>1.d-3) write(19, '("Qs[", i3, "]=", f13.3)') i, xmin(i, 1)     
        !     enddo
        !     call cpu_time(finish)
        !     print *, log_filename
        !     j = 0
        !     do i=1, n_hex
        !         if(abs(xmin(i, 1))>1.d-3) then
        !             j = j + 1
        !             write(20, '("Qs[", i3, "]=", f13.3)') i, xmin(i, 1)     
        !         endif
        !     enddo
        !     n_hex_global = j
        !     y_min_global = ymin
        !     write(*, '("ymin =", f15.3)') ymin
        !     write(*, '("Time =", f10.1)') finish-start
        !     write(20, '("ymin =", f15.3)') ymin
        !     write(20, '("Time =", f10.1)') finish-start
        !     close(20)
 
        !     return
        ! end subroutine evolve_fixed

        ! subroutine run_fixed_sade(n, idx_q, case_name, stage, np, max_iter, &
        !         random_state, qmin, lp, prob)
        !     implicit none
        !     character(len=*) :: case_name
        !     integer(kind=4) :: n, idx_q(n), stage, max_iter, np, lp
        !     real(kind=8) :: qmin, random_state, prob
        !     integer(kind=4) :: i, j
        !     logical :: exist
        !     real(kind=8) :: r

        !     call init_de(case_name, stage, np)
        !     call random_number(r)
        !     print *, r
        !     inquire(file="log/run_fix_sade.txt", exist=exist)
        !     if (exist) then
        !       open(13, file="log/run_fix_sade.txt", status="old", position="append", action="write")
        !     else
        !       open(13, file="log/run_fix_sade.txt", status="new", action="write")
        !     end if
        !     rn(1) = r
        !     call init_fixed_population(n, idx_q, np)
        !     call evolve_fixed(n, idx_q, np, max_iter, qmin, lp, prob)
        !     write(13, *) log_filename, r, n_hex_global, y_min_global 
        !     close(13)

        ! end subroutine run_fixed_sade

end module sade 
