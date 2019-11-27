program sade_9sp1
    use mpi_sade
    implicit none
    call mpi_run_sade("9sp1", 5, 5*5*4*10, 500000, 100.d0, 100, 100)
end program sade_9sp1
