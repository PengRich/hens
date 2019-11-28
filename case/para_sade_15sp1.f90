program sade_15sp1
    use mpi_sade
    implicit none
    call run_parallel_sade("15sp1", 2, 2*7*8*10, 10000, 250.d0, 100, 100, 30, &
        0.3d0, 0.3d0)
end program sade_15sp1
