program sade_9sp1
    use mpi_sade
    implicit none
    call run_parallel_sade("4sp1", 2, 2*2*2*10, 100, 100.d0, 100, 100, 5, 10, 0.025d0)
end program sade_9sp1
