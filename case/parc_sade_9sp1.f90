program sade_9sp1
    use mpi_sade
    implicit none
    call run_parallel_sade_c("9sp1", 4, 4*5*4*10, 10000, 250.d0, 100, 100, 30)
end program sade_9sp1
