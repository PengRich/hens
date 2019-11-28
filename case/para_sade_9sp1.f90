program sade_9sp1
    use mpi_sade
    implicit none
    call run_parallel_sade("9sp1", 4, 4*5*4*10, 10000, 250.d0, 100, 100, 30, &
        0.3d0, 0.3d0)
end program sade_9sp1
