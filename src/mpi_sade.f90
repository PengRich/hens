module mpi_sade
    use mpi
    use sade
    implicit none
    contains
        subroutine mpi_run_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number
            real(kind=8), intent(in) :: qmin !, elimination_prob

            logical :: exist
            character(len=23) :: filename
            real(kind=8) :: random_state
            ! MPI
            integer :: node, n_core, ierr, status(mpi_status_size), i
            logical :: Ionode
            character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
            integer :: namelen
            character(len=1) :: node0
            character(len=10) :: perfix

            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, n_core, ierr)
            call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)

            call cpu_time(start)
            n_open = node + 20
            n_log_file = n_open + 20

            write(node0, '(i1)') node

            filename = "output/"// node0 // "_" // trim(case_name) // "_sade.txt"
            inquire(file=filename, exist=exist)
            if (exist) then
                open(n_open, file=filename, status="old", position="append", action="write")
            else
                open(n_open, file=filename, status="new", action="write")
            end if

            Ionode=(node .eq. 0)
            if(Ionode) then
                call set_random_seed()
                do i=1, n_core-1
                    do while(.true.)
                        call random_number(random_state)
                        if(random_state > 0.1) exit
                    enddo
                    call MPI_SSEND(random_state, 1, MPI_DOUBLE_PRECISION, i, 99, &
                            MPI_COMM_WORLD, ierr)
                enddo
                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1) exit
                enddo
            else
                call MPI_RECV(random_state, 1, MPI_DOUBLE_PRECISION, 0, 99, &
                    MPI_COMM_WORLD, status, ierr)
            endif
            rn(1) = dble(random_state)
            perfix = trim(node0 // "_" // case_name)
            call get_log_filename(perfix)
            open(unit=n_log_file, file=log_filename, action="write", status="replace")
            call init_de(case_name, stage, np)
            call init_population(np)
            call init_de_parameter(np, sampling_number)
            call evolve(np, max_iter, qmin, learning_period)
            call deallocate_de_parameter()
            call deallocate_de_var()
            call deallocate_var()
            call cpu_time(finish)
            write(*, *) log_filename, dble(random_state), n_hex_global, y_min_global
            write(n_open, *) log_filename, dble(random_state), n_hex_global, &
                y_min_global, finish-start
            close(n_open)
            close(n_log_file)
            call MPI_FINALIZE(ierr)
            return
        end subroutine mpi_run_sade

        subroutine run_parallel_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, switch_number, switch_ratio, &
                reinit_ratio)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number, switch_number
            real(kind=8), intent(in) :: qmin, switch_ratio, reinit_ratio

            logical :: exist
            character(len=28) :: filename
            ! MPI
            integer :: node, n_core, ierr, status(mpi_status_size), i, &
                node_trg, node_src
            logical :: Ionode
            character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
            integer :: namelen
            character(len=1) :: node0
            character(len=15) :: perfix
            ! Parallel
            real(kind=8) :: random_state
            integer(kind=4) :: j, k, idx_np, idx_best(1), &
                iter_switch, n_switch, n_recv

            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, n_core, ierr)
            call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)

            call cpu_time(start)
            print *, "start process", node
            write(node0, '(i1)') node

            filename = "output/"// node0 // "_para_" // trim(case_name) // "_sade.txt"
            inquire(file=filename, exist=exist)
            if (exist) then
                open(n_open, file=filename, status="old", position="append", action="write")
            else
                open(n_open, file=filename, status="new", action="write")
            end if
            close(n_open)

            n_open = node + 20
            n_log_file = n_open + 20
            iter_switch = switch_number*n_core

            Ionode=(node .eq. 0)
            if(Ionode) then
                call set_random_seed()
                do i=1, n_core-1
                    do while(.true.)
                        call random_number(random_state)
                        if(random_state > 0.1) exit
                    enddo
                    call MPI_SSEND(random_state, 1, MPI_DOUBLE_PRECISION, i, 99, &
                            MPI_COMM_WORLD, ierr)
                enddo
                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1) exit
                enddo
            else
                call MPI_RECV(random_state, 1, MPI_DOUBLE_PRECISION, 0, 99, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            rn(1) = dble(random_state)

            call init_de(case_name, stage, np)
            call init_population(np)
            call init_de_parameter(np, sampling_number)

            n_switch = int(switch_ratio * real(np))
            perfix = trim(node0 // "_para_" // case_name)
            do j=1, iter_switch
                open(n_open, file=filename, status="old", position="append", action="write")
                call get_log_filename(perfix)
                open(unit=n_log_file, file=log_filename, action="write", status="replace")
                call evolve(np, max_iter, qmin, learning_period)
                call cpu_time(finish)
                close(n_log_file)
                write(*, *) iter_switch-j+1, node, log_filename, n_hex_global, &
                    y_min_global, finish-start
                write(n_open, *) iter_switch-j+1, log_filename, n_hex_global, &
                    y_min_global, finish-start
                close(n_open)

                if(node<n_core-1) then
                    node_trg = node + 1
                else
                    node_trg = 0
                endif
                if(node==0) then
                    node_src = n_core - 1
                else
                    node_src = node - 1
                endif

                idx_best = minloc(y_old)

                call MPI_SEND(x(:, idx_best(1)), n_hex, MPI_DOUBLE_PRECISION, &
                    node_trg, 99, MPI_COMM_WORLD, ierr)
                call MPI_RECV(x(:, idx_best(1)), n_hex, MPI_DOUBLE_PRECISION, &
                    node_src, 99, MPI_COMM_WORLD, status, ierr)
                y_old(idx_best(1)) = tac(x(:, idx_best(1)))
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                n_recv = 0
                do k=1, n_switch
                    do while(.true.)
                        idx_np = int(rand(rn(1))*np)+1
                        if(idx_np==idx_best(1)) cycle
                        exit
                    enddo
                    ! idx_np = int(rand(rn(1))*np)+1
                    call MPI_SEND(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_trg, 99, MPI_COMM_WORLD, ierr)
                    ! idx_np = int(rand(rn(1))*np)+1
                    call MPI_RECV(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_src, 99, MPI_COMM_WORLD, status, ierr)
                    y_old(idx_np) = tac(x(:, idx_np))
                    do i=1, n_hex
                        if(x(i, idx_np)>0.d0) n_recv = n_recv + 1
                    enddo
                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                enddo
                call reinit_population(np, idx_best(1), reinit_ratio)
                ! ymin = minval(y_old)
                ! xmin = x(:, minloc(y_old))
                ! print *, "node", node, "completed"
                print *, node, n_recv, n_switch, n_recv/n_switch
            enddo

            call deallocate_de_parameter()
            call deallocate_de_var()
            call deallocate_var()

            call cpu_time(finish)
            open(n_open, file=filename, status="old", position="append", action="write")
            write(n_open, *) dble(random_state), finish-start
            close(n_open)

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if(Ionode) then
                filename = "output/"// "s_para_" // trim(case_name) // "_sade.txt"
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(n_open, file=filename, status="old", position="append", action="write")
                else
                    open(n_open, file=filename, status="new", action="write")
                end if
                write(n_open, *) 0, y_min_global
                do i=1, n_core-1
                    call MPI_RECV(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                        i, 99, MPI_COMM_WORLD, status, ierr)
                    write(n_open, *) i, y_min_global
                enddo
                close(n_open)
            else
                call MPI_SEND(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                    0, 99, MPI_COMM_WORLD, ierr)
            endif
 
            call MPI_FINALIZE(ierr)
            return
        end subroutine run_parallel_sade

       subroutine run_parallel_sade_b(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, switch_number, switch_ratio, &
                reinit_ratio)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number, switch_number
            real(kind=8), intent(in) :: qmin, switch_ratio, reinit_ratio

            logical :: exist
            character(len=28) :: filename
            ! MPI
            integer :: node, n_core, ierr, status(mpi_status_size), i, &
                node_trg, node_src
            logical :: Ionode
            character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
            integer :: namelen
            character(len=1) :: node0
            character(len=15) :: perfix
            ! Parallel
            real(kind=8) :: random_state
            integer(kind=4) :: j, k, idx_np, idx_best(1), &
                iter_switch, n_switch, n_recv

            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, n_core, ierr)
            call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)

            call cpu_time(start)
            print *, "start process", node
            write(node0, '(i1)') node

            filename = "output/"// node0 // "_para_" // trim(case_name) // "_sade.txt"
            inquire(file=filename, exist=exist)
            if (exist) then
                open(n_open, file=filename, status="old", position="append", action="write")
            else
                open(n_open, file=filename, status="new", action="write")
            end if
            close(n_open)

            n_open = node + 20
            n_log_file = n_open + 20
            iter_switch = switch_number*n_core

            Ionode=(node .eq. 0)
            if(Ionode) then
                call set_random_seed()
                do i=1, n_core-1
                    do while(.true.)
                        call random_number(random_state)
                        if(random_state > 0.1) exit
                    enddo
                    call MPI_SSEND(random_state, 1, MPI_DOUBLE_PRECISION, i, 99, &
                            MPI_COMM_WORLD, ierr)
                enddo
                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1) exit
                enddo
            else
                call MPI_RECV(random_state, 1, MPI_DOUBLE_PRECISION, 0, 99, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            rn(1) = dble(random_state)

            call init_de(case_name, stage, np)
            call init_population(np)
            call init_de_parameter(np, sampling_number)

            n_switch = int(switch_ratio * real(np))
            perfix = trim(node0 // "_para_" // case_name)
            do j=1, iter_switch
                open(n_open, file=filename, status="old", position="append", action="write")
                call get_log_filename(perfix)
                open(unit=n_log_file, file=log_filename, action="write", status="replace")
                call evolve(np, max_iter, qmin, learning_period)
                close(n_log_file)
                write(*, *) node, log_filename, n_hex_global, y_min_global
                write(n_open, *) log_filename, n_hex_global, y_min_global
                close(n_open)

                if(node<n_core-1) then
                    node_trg = node + 1
                else
                    node_trg = 0
                endif
                if(node==0) then
                    node_src = n_core - 1
                else
                    node_src = node - 1
                endif

                n_recv = 0
                do k=1, n_switch
                    idx_np = int(rand(rn(1))*np)+1
                    call MPI_SEND(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_trg, 99, MPI_COMM_WORLD, ierr)
                    call MPI_RECV(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_src, 99, MPI_COMM_WORLD, status, ierr)
                    y_old(idx_np) = tac(x(:, idx_np))
                    do i=1, n_hex
                        if(x(i, idx_np)>0.d0) n_recv = n_recv + 1
                    enddo
                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                enddo

                idx_best = minloc(y_old)
                call reinit_population(np, idx_best(1), reinit_ratio)
                ! ymin = minval(y_old)
                ! xmin = x(:, minloc(y_old))
                print *, node, n_recv, n_switch, n_recv/n_switch
            enddo

            call deallocate_de_parameter()
            call deallocate_de_var()
            call deallocate_var()

            call cpu_time(finish)
            open(n_open, file=filename, status="old", position="append", action="write")
            write(n_open, *) dble(random_state), finish-start
            close(n_open)

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if(Ionode) then
                filename = "output/"// "s_para_" // trim(case_name) // "_sade.txt"
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(n_open, file=filename, status="old", position="append", action="write")
                else
                    open(n_open, file=filename, status="new", action="write")
                end if
                write(n_open, *) 0, y_min_global
                do i=1, n_core-1
                    call MPI_RECV(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                        i, 99, MPI_COMM_WORLD, status, ierr)
                    write(n_open, *) i, y_min_global
                enddo
                close(n_open)
            else
                call MPI_SEND(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                    0, 99, MPI_COMM_WORLD, ierr)
            endif

            call MPI_FINALIZE(ierr)
            return
        end subroutine run_parallel_sade_b

        subroutine run_parallel_sade_c(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, iter_switch)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number
            real(kind=8), intent(in) :: qmin

            logical :: exist
            character(len=28) :: filename
            ! MPI
            integer :: node, n_core, ierr, status(mpi_status_size), i, &
                node_trg, node_src
            logical :: Ionode
            character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
            integer :: namelen
            character(len=1) :: node0
            character(len=15) :: perfix
            ! Parallel
            real(kind=8) :: random_state
            integer(kind=4) :: j, k, idx_np, idx_best(1), idx_worst(1), iter_switch

            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, n_core, ierr)
            call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)

            call cpu_time(start)
            print *, "start process", node
            write(node0, '(i1)') node

            filename = "output/"// node0 // "_para_" // trim(case_name) // "_sade.txt"
            inquire(file=filename, exist=exist)
            if (exist) then
                open(n_open, file=filename, status="old", position="append", action="write")
            else
                open(n_open, file=filename, status="new", action="write")
            end if
            close(n_open)

            n_open = node + 20
            n_log_file = n_open + 20

            Ionode=(node .eq. 0)
            if(Ionode) then
                call set_random_seed()
                do i=1, n_core-1
                    do while(.true.)
                        call random_number(random_state)
                        if(random_state > 0.1) exit
                    enddo
                    call MPI_SSEND(random_state, 1, MPI_DOUBLE_PRECISION, i, 99, &
                            MPI_COMM_WORLD, ierr)
                enddo
                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1) exit
                enddo
            else
                call MPI_RECV(random_state, 1, MPI_DOUBLE_PRECISION, 0, 99, &
                    MPI_COMM_WORLD,status,ierr)
            endif
            rn(1) = dble(random_state)

            call init_de(case_name, stage, np)
            call init_population(np)
            call init_de_parameter(np, sampling_number)

            perfix = trim(node0 // "_parc_" // case_name)
            do j=1, iter_switch
                open(n_open, file=filename, status="old", position="append", action="write")
                call get_log_filename(perfix)
                open(unit=n_log_file, file=log_filename, action="write", status="replace")
                call evolve(np, max_iter, qmin, learning_period)
                close(n_log_file)
                write(*, *) node, log_filename, n_hex_global, y_min_global
                write(n_open, *) log_filename, n_hex_global, y_min_global
                close(n_open)

                if(node<n_core-1) then
                    node_trg = node + 1
                else
                    node_trg = 0
                endif
                if(node==0) then
                    node_src = n_core - 1
                else
                    node_src = node - 1
                endif

                if(Ionode) then
                    do k=1, n_core-1
                        idx_worst = maxloc(y_old)
                        print *, k, "before", tac(x(:, idx_worst(1)))
                        call MPI_RECV(x(:, idx_worst(1)), n_hex, MPI_DOUBLE_PRECISION, &
                            k, 99, MPI_COMM_WORLD, status, ierr)
                        print *, k, "after", tac(x(:, idx_worst(1)))
                        y_old(idx_worst(1)) = tac(x(:, idx_worst(1)))
                    enddo
                    ymin = minval(y_old)
                    print *, "global", ymin
                    xmin = x(:, minloc(y_old))
                else
                    idx_best = minloc(y_old)
                    call MPI_SEND(x(:, idx_best(1)), n_hex, MPI_DOUBLE_PRECISION, &
                        0, 99, MPI_COMM_WORLD, ierr)
                    call reinit_population(np, 0, 1.d0)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            enddo

            call deallocate_de_parameter()
            call deallocate_de_var()
            call deallocate_var()

            call cpu_time(finish)
            open(n_open, file=filename, status="old", position="append", action="write")
            write(n_open, *) dble(random_state), finish-start
            close(n_open)

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if(Ionode) then
                filename = "output/"// "s_parc_" // trim(case_name) // "_sade.txt"
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(n_open, file=filename, status="old", position="append", action="write")
                else
                    open(n_open, file=filename, status="new", action="write")
                end if
                write(n_open, *) 0, y_min_global
                do i=1, n_core-1
                    call MPI_RECV(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                        i, 99, MPI_COMM_WORLD, status, ierr)
                    write(n_open, *) i, y_min_global
                enddo
                close(n_open)
            else
                call MPI_SEND(y_min_global, 1, MPI_DOUBLE_PRECISION, &
                    0, 99, MPI_COMM_WORLD, ierr)
            endif

            call MPI_FINALIZE(ierr)
            return
        end subroutine run_parallel_sade_c

end module mpi_sade
