program sade_9sp1
    use mpi_sade
    implicit none
    call run_parallel_sade("10sp1", 4, 5*5*4*10, 100, 100.d0, 100, 100, 3, 0.d0)
end program sade_9sp1
