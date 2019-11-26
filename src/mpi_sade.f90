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
            ! MPI
            real :: random_state
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
                    call MPI_SSEND(random_state, 1, MPI_INTEGER, i, 99, &
                            MPI_COMM_WORLD, ierr)
                enddo
                do while(.true.)
                    call random_number(random_state)
                    if(random_state > 0.1) exit
                enddo
            else
                call MPI_RECV(random_state, 1, MPI_INTEGER, 0, 99, &
                    MPI_COMM_WORLD,status,ierr)
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

end module mpi_sade
