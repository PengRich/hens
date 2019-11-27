module mpi_sade
    use mpi
    use sade
    implicit none
    contains
        subroutine mpi_run_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, elimination_number)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number, elimination_number
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
            call evolve(np, max_iter, qmin, learning_period, &
                sampling_number, elimination_number)
            call deallocate_de_var()
            call deallocate_var()
            write(*, *) log_filename, dble(random_state), n_hex_global, y_min_global
            write(n_open, *) log_filename, dble(random_state), n_hex_global, y_min_global
            close(n_open)
            close(n_log_file)
            call MPI_FINALIZE(ierr)
            return
        end subroutine mpi_run_sade

        subroutine run_parallel_sade(case_name, stage, np, max_iter, qmin, &
                learning_period, sampling_number, elimination_number, &
                switch_number, switch_ratio)
            implicit none
            character(len=*), intent(in) :: case_name
            integer(kind=4), intent(in) :: stage, max_iter, np, &
                learning_period, sampling_number, elimination_number, &
                switch_number
            real(kind=8), intent(in) :: qmin, switch_ratio !, elimination_prob

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
            integer(kind=4) :: j, k, idx_np, iter_switch, n_switch

            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, node, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, n_core, ierr)
            call MPI_GET_PROCESSOR_NAME(hostname, namelen, ierr)
            print *, "start process", node
            n_open = node + 20
            n_log_file = n_open + 20
            iter_switch = switch_number

            write(node0, '(i1)') node

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

            n_switch = int(switch_ratio * real(np))
            do j=1, iter_switch
                filename = "output/"// node0 // "_para_" // trim(case_name) // "_sade.txt"
                inquire(file=filename, exist=exist)
                if (exist) then
                    open(n_open, file=filename, status="old", position="append", action="write")
                else
                    open(n_open, file=filename, status="new", action="write")
                end if

                perfix = trim(node0 // "_para_" // case_name)
                call get_log_filename(perfix)
                open(unit=n_log_file, file=log_filename, action="write", status="replace")
                call evolve(np, max_iter, qmin, learning_period, &
                    sampling_number, elimination_number)
                close(n_log_file)

                write(*, *) node, log_filename, dble(random_state), n_hex_global, y_min_global
                write(n_open, *) log_filename, dble(random_state), n_hex_global, y_min_global
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

                call MPI_SEND(xmin(:, 1), n_hex, MPI_DOUBLE_PRECISION, &
                    node_trg, 99, MPI_COMM_WORLD, ierr)
                idx_np = int(rand(rn(1))*np)+1
                call MPI_RECV(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                    node_src, 99, MPI_COMM_WORLD, status, ierr)

                do k=1, n_switch
                    idx_np = int(rand(rn(1))*np)+1
                    call MPI_SEND(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_trg, 99, MPI_COMM_WORLD, ierr)
                    idx_np = int(rand(rn(1))*np)+1
                    call MPI_RECV(x(:, idx_np), n_hex, MPI_DOUBLE_PRECISION, &
                        node_src, 99, MPI_COMM_WORLD, status, ierr)
                enddo

                ymin = minval(y_old)
                xmin = x(:, minloc(y_old))

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            enddo

            call deallocate_de_var()
            call deallocate_var()
            print *, "node", node, "completed"
            call MPI_FINALIZE(ierr)
            return
        end subroutine run_parallel_sade

end module mpi_sade
