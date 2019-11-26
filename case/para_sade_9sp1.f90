program sade_9sp1
    use mpi_sade
    implicit none
    call run_parallel_sade("9sp1", 5, 5*5*4*10, 100000, 100.d0, 100, 100, 15, 5, 0.d0)
end program sade_9sp1
